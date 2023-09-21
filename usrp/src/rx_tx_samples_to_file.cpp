//
// Copyright 2010-2012,2014-2015 Ettus Research LLC
// Copyright 2018 Ettus Research, a National Instruments Company
//
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <csignal>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>

namespace po = boost::program_options;

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
static bool start_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

/***********************************************************************
 * Utilities
 **********************************************************************/
//! Change to filename, e.g. from usrp_samples.dat to usrp_samples.00.dat,
//  but only if multiple names are to be generated.
std::string generate_out_filename(
    const std::string& base_fn, size_t n_names, size_t this_name)
{
    if (n_names == 1) {
        return base_fn;
    }

    boost::filesystem::path base_fn_fp(base_fn);
    base_fn_fp.replace_extension(boost::filesystem::path(
        str(boost::format("%02d%s") % this_name % base_fn_fp.extension().string())));
    return base_fn_fp.string();
}


/***********************************************************************
 * transmit_worker function
 * A function to be used in a thread for transmitting
 **********************************************************************/
template <typename samp_type>
void send_from_file(
    uhd::usrp::multi_usrp::sptr usrp, 
    uhd::tx_streamer::sptr tx_stream, 
    const std::string& file, 
    size_t samps_per_buff, 
    uhd::time_spec_t time_spec)
{
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst   = false;
    std::vector<samp_type> buff(samps_per_buff,0);
    std::ifstream infile(file.c_str(), std::ifstream::binary);

    tx_stream->send("", 0, md);

    // loop until the entire file has been read
    while (not md.end_of_burst and not stop_signal_called) {
        infile.read((char*)&buff.front(), buff.size() * sizeof(samp_type));
        size_t num_tx_samps = size_t(infile.gcount() / sizeof(samp_type));

        //std::cout << "Num tx samps: " << num_tx_samps << std::endl;

        md.end_of_burst = infile.eof();

        if(md.start_of_burst) {
            md.has_time_spec = true;
            md.time_spec = time_spec;
            std::cout << "Set settling time" << std::endl;

            //usrp->set_command_time(uhd::time_spec_t(settling_time));
        }
        else {
            md.has_time_spec = false;
        }

       // while(not start_signal_called) {} // wait for start signal
        for (size_t sent_all = 0;	 sent_all != num_tx_samps; ) {
		    sent_all += tx_stream->send(&buff.front() + sent_all, num_tx_samps - sent_all, md);
        }


        md.start_of_burst = false;
    }

    // send a mini EOB packet
    std::fill(buff.begin(),buff.end(),0);
    md.end_of_burst = true;
    tx_stream->send(&buff.front(), samps_per_buff, md);

    // added to allow tx streamer to still read buffer before closing file
    //std::this_thread::sleep_for(std::chrono::seconds(1));

    std::cout << "TX done!" << std::endl;

    infile.close();
}

typedef std::function<uhd::sensor_value_t(const std::string&)> get_sensor_fn_t;

// function that blocks until USRP is seen to be locked to reference
bool check_locked_sensor(std::vector<std::string> sensor_names,
    const char* sensor_name,
    get_sensor_fn_t get_sensor_fn,
    double setup_time)
{
    if (std::find(sensor_names.begin(), sensor_names.end(), sensor_name)
        == sensor_names.end())
        return false;

    auto setup_timeout = std::chrono::steady_clock::now()
                         + std::chrono::milliseconds(int64_t(setup_time * 1000));
    bool lock_detected = false;

    std::cout << boost::format("Waiting for \"%s\": ") % sensor_name;
    std::cout.flush();

    while (true) {
        if (lock_detected) {
            std::cout << " locked." << std::endl;
            break;
        }
        if (get_sensor_fn(sensor_name).to_bool()) {
            std::cout << "+";
            std::cout.flush();
            lock_detected = true;
        } else {
            if (std::chrono::steady_clock::now() > setup_timeout) {
                std::cout << std::endl;
                throw std::runtime_error(
                    str(boost::format(
                            "timed out waiting for consecutive locks on sensor \"%s\"")
                        % sensor_name));
            }
            std::cout << "_";
            std::cout.flush();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    std::cout << std::endl;
    return true;
}

/***********************************************************************
 * recv_to_file function
 **********************************************************************/
template <typename samp_type>
void recv_to_file(uhd::usrp::multi_usrp::sptr usrp,
    const std::string& cpu_format,
    const std::string& wire_format,
    const std::string& file,
    size_t samps_per_buff,
    int num_requested_samples,
    std::vector<size_t> rx_channel_nums,
    double settling_time,
    uhd::time_spec_t time_spec)
{
    int num_total_samps = 0;
    // create a receive streamer
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    stream_args.channels             = rx_channel_nums;
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    // Prepare buffers for received samples and metadata
    uhd::rx_metadata_t md;
    std::vector<std::vector<samp_type>> buffs(
        rx_channel_nums.size(), std::vector<samp_type>(samps_per_buff));
    // create a vector of pointers to point to each of the channel buffers
    std::vector<samp_type*> buff_ptrs;
    for (size_t i = 0; i < buffs.size(); i++) {
        for (size_t x = 0; x < samps_per_buff; x++) {
            buffs[i][x] = 0;
        }
        buff_ptrs.push_back(&buffs[i].front());
    }

    // Create one ofstream object per channel
    // (use shared_ptr because ofstream is non-copyable)
    std::vector<std::shared_ptr<std::ofstream>> outfiles;
    for (size_t i = 0; i < buffs.size(); i++) {
        const std::string this_filename = generate_out_filename(file, buffs.size(), i);
        outfiles.push_back(std::shared_ptr<std::ofstream>(
            new std::ofstream(this_filename.c_str(), std::ofstream::binary)));
    }
    UHD_ASSERT_THROW(outfiles.size() == buffs.size());
    UHD_ASSERT_THROW(buffs.size() == rx_channel_nums.size());
    bool overflow_message = true;
    // We increase the first timeout to cover for the delay between now + the
    // command time, plus 500ms of buffer. In the loop, we will then reduce the
    // timeout for subsequent receives.
    double timeout = settling_time + 0.5f;

    // setup streaming
    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)
                                     ? uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS
                                     : uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    stream_cmd.num_samps  = num_requested_samples;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec  = time_spec;
    rx_stream->issue_stream_cmd(stream_cmd);
    //start_signal_called = true;

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)) {
        size_t num_rx_samps = rx_stream->recv(buff_ptrs, samps_per_buff, md, timeout);
        timeout             = 0.1f; // small timeout for subsequent recv

        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << "Timeout while streaming" << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
            if (overflow_message) {
                overflow_message = false;
                std::cerr
                    << boost::format(
                           "Got an overflow indication. Please consider the following:\n"
                           "  Your write medium must sustain a rate of %fMB/s.\n"
                           "  Dropped samples will not be written to the file.\n"
                           "  Please modify this example for your purposes.\n"
                           "  This message will not appear again.\n")
                           % (usrp->get_rx_rate() * sizeof(samp_type) / 1e6);
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
            throw std::runtime_error("Receiver error " + md.strerror());
        }

        num_total_samps += num_rx_samps;

        for (size_t i = 0; i < outfiles.size(); i++) {
            outfiles[i]->write(
                (const char*)buff_ptrs[i], num_rx_samps * sizeof(samp_type));
        }
    }

    // Shut down receiver
    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);

    // Close files
    for (size_t i = 0; i < outfiles.size(); i++) {
        outfiles[i]->close();
    }
}


/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char* argv[])
{
    // transmit variables to be set by po
    std::string args, tx_subdev, ref, otw, tx_channels;
    double tx_rate;

    // receive variables to be set by po
    std::string rx_file, tx_file, type, rx_subdev, rx_channels, sync;
    size_t total_num_samps, spb;
    double rx_rate, settling, lock_tm;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "uhd device address args")
        //("rx-args", po::value<std::string>(&rx_args)->default_value(""), "uhd receive device address args")
        ("tx-file", po::value<std::string>(&tx_file)->default_value("tx_usrp_samples.dat"), "name of the file to read binary samples from")
        ("rx-file", po::value<std::string>(&rx_file)->default_value("rx_usrp_samples.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("float"), "sample type in file: double, float, or short")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("settling", po::value<double>(&settling)->default_value(double(0.2)), "settling time (seconds) before receiving")
        ("spb", po::value<size_t>(&spb)->default_value(0), "samples per buffer, 0 for default")
        ("tx-rate", po::value<double>(&tx_rate), "rate of transmit outgoing samples")
        ("rx-rate", po::value<double>(&rx_rate), "rate of receive incoming samples")
        // ("tx-freq", po::value<double>(&tx_freq), "transmit RF center frequency in Hz")
        // ("rx-freq", po::value<double>(&rx_freq), "receive RF center frequency in Hz")
        // ("ampl", po::value<float>(&ampl)->default_value(float(0.3)), "amplitude of the waveform [0 to 0.7]")
        // ("tx-gain", po::value<double>(&tx_gain), "gain for the transmit RF chain")
        // ("rx-gain", po::value<double>(&rx_gain), "gain for the receive RF chain")
        // ("tx-ant", po::value<std::string>(&tx_ant), "transmit antenna selection")
        // ("rx-ant", po::value<std::string>(&rx_ant), "receive antenna selection")
        ("tx-subdev", po::value<std::string>(&tx_subdev)->default_value("A:AB"), "transmit subdevice specification")
        ("rx-subdev", po::value<std::string>(&rx_subdev)->default_value("A:AB"), "receive subdevice specification")
        // ("tx-bw", po::value<double>(&tx_bw), "analog transmit filter bandwidth in Hz")
        // ("rx-bw", po::value<double>(&rx_bw), "analog receive filter bandwidth in Hz")
        // ("wave-type", po::value<std::string>(&wave_type)->default_value("CONST"), "waveform type (CONST, SQUARE, RAMP, SINE)")
        // ("wave-freq", po::value<double>(&wave_freq)->default_value(0), "waveform frequency in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("external"), "clock reference (internal, external, mimo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode")
        ("tx-channels", po::value<std::string>(&tx_channels)->default_value("0"), "which TX channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("rx-channels", po::value<std::string>(&rx_channels)->default_value("0"), "which RX channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("sync", po::value<std::string>(&sync)->default_value("pps"), "time sync method. Specify pps or now")
        ("lock-tm", po::value<double>(&lock_tm)->default_value(0.1), "time to wait in seconds before lock is detected, recommend 0.1")
        // ("tx-int-n", "tune USRP TX with integer-N tuning")
        // ("rx-int-n", "tune USRP RX with integer-N tuning")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << "UHD TXRX Samples to File " << desc << std::endl;
        return ~0;
    }

    // create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the transmit/receive usrp device with: %s...") % args
              << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
    // std::cout << std::endl;
    // std::cout << boost::format("Creating the receive usrp device with: %s...") % rx_args
    //           << std::endl;
    // uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

    // detect which channels to use
    std::vector<std::string> tx_channel_strings;
    std::vector<size_t> tx_channel_nums;
    boost::split(tx_channel_strings, tx_channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < tx_channel_strings.size(); ch++) {
        size_t chan = std::stoi(tx_channel_strings[ch]);
        if (chan >= usrp->get_tx_num_channels()) {
            throw std::runtime_error("Invalid TX channel(s) specified.");
        } else
            tx_channel_nums.push_back(std::stoi(tx_channel_strings[ch]));
    }
    std::vector<std::string> rx_channel_strings;
    std::vector<size_t> rx_channel_nums;
    boost::split(rx_channel_strings, rx_channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < rx_channel_strings.size(); ch++) {
        size_t chan = std::stoi(rx_channel_strings[ch]);
        if (chan >= usrp->get_rx_num_channels()) {
            throw std::runtime_error("Invalid RX channel(s) specified.");
        } else
            rx_channel_nums.push_back(std::stoi(rx_channel_strings[ch]));
    }

    // always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("tx-subdev"))
        usrp->set_tx_subdev_spec(tx_subdev,tx_channel_nums[0]);
    if (vm.count("rx-subdev"))
        usrp->set_rx_subdev_spec(rx_subdev,rx_channel_nums[0]);

    // Lock mboard clocks
    if (ref == "external") {
        usrp->set_clock_source(ref);
        // rx_usrp->set_clock_source(ref);
    }

    std::cout << "Using Devices: " << usrp->get_pp_string() << std::endl;
    //std::cout << "Using RX Device: " << rx_usrp->get_pp_string() << std::endl;

    // set the transmit sample rate
    if (not vm.count("tx-rate")) {
        std::cerr << "Please specify the transmit sample rate with --tx-rate"
                  << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate / 1e6)
              << std::endl;
    usrp->set_tx_rate(tx_rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...")
                     % (usrp->get_tx_rate() / 1e6)
              << std::endl
              << std::endl;

    // set the receive sample rate
    if (not vm.count("rx-rate")) {
        std::cerr << "Please specify the sample rate with --rx-rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rx_rate / 1e6)
              << std::endl;
    usrp->set_rx_rate(rx_rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...")
                     % (usrp->get_rx_rate() / 1e6)
              << std::endl
              << std::endl;

    // create a transmit streamer
    // linearly map channels (index0 = channel0, index1 = channel1, ...)
    uhd::stream_args_t stream_args("fc32", otw);
    stream_args.channels             = tx_channel_nums;
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

    // allocate a buffer which we re-use for each channel
    if (spb == 0)
        spb = tx_stream->get_max_num_samps();
    //sstd::vector<std::complex<float>> buff(spb);
    // int num_channels = tx_channel_nums.size();

    // // Check Ref and LO Lock detect
    
    // tx_sensor_names = usrp->get_tx_sensor_names(0);
    // if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked")
    //     != tx_sensor_names.end()) {
    //     uhd::sensor_value_t lo_locked = usrp->get_tx_sensor("lo_locked", 0);
    //     std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string()
    //               << std::endl;
    //     UHD_ASSERT_THROW(lo_locked.to_bool());
    // }
    // rx_sensor_names = usrp->get_rx_sensor_names(0);
    // if (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "lo_locked")
    //     != rx_sensor_names.end()) {
    //     uhd::sensor_value_t lo_locked = usrp->get_rx_sensor("lo_locked", 0);
    //     std::cout << boost::format("Checking RX: %s ...") % lo_locked.to_pp_string()
    //               << std::endl;
    //     UHD_ASSERT_THROW(lo_locked.to_bool());
    // }
    
    // std::vector<std::string> tx_sensor_names, rx_sensor_names;
    // tx_sensor_names = usrp->get_mboard_sensor_names(0);
    // rx_sensor_names = usrp->get_mboard_sensor_names(0);


    // if ((ref == "mimo")
    //     and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "mimo_locked")
    //             != tx_sensor_names.end())) {
    //     uhd::sensor_value_t mimo_locked = usrp->get_mboard_sensor("mimo_locked", 0);
    //     std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string()
    //               << std::endl;
    //     UHD_ASSERT_THROW(mimo_locked.to_bool());
    // }
    size_t num_mboards = usrp->get_num_mboards();

    if (ref == "external") {

        // check RX reference
        for(size_t mboard = 0; mboard < num_mboards; ++mboard) {
            std::cout << boost::format("Locking mboard %d...") % mboard << std::endl;
            // uhd::sensor_value_t ref_locked = usrp->get_mboard_sensor("ref_locked", mboard);
            // std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string()
            //       << std::endl;
            // UHD_ASSERT_THROW(ref_locked.to_bool());

            check_locked_sensor(usrp->get_mboard_sensor_names(mboard),
                "ref_locked",
                [usrp, mboard](const std::string& sensor_name) {
                    return usrp->get_mboard_sensor(sensor_name, mboard);
                },lock_tm);
        }
    }

    // if ((ref == "mimo")
    //     and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "mimo_locked")
    //             != rx_sensor_names.end())) {
    //     uhd::sensor_value_t mimo_locked = usrp->get_mboard_sensor("mimo_locked", 0);
    //     std::cout << boost::format("Checking RX: %s ...") % mimo_locked.to_pp_string()
    //               << std::endl;
    //     UHD_ASSERT_THROW(mimo_locked.to_bool());
    // }

        // uhd::sensor_value_t ref_locked = usrp->get_mboard_sensor("ref_locked", 0);
        // std::cout << boost::format("Checking RX: %s ...") % ref_locked.to_pp_string()
        //           << std::endl;
        // UHD_ASSERT_THROW(ref_locked.to_bool());

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    // reset usrp time to prepare for transmit/receive
    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;

    // set pps
    if(sync == "pps") {
        usrp->set_time_source("external");
        usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
        std::this_thread::sleep_for(std::chrono::seconds(1)); // wait for pps sync pulse
    }
    else if(sync == "now") {
        // This is not a true time lock, the devices will be off by a few RTT.
        // Rather, this is just to allow for demonstration of the code below.
        usrp->set_time_now(uhd::time_spec_t(0.0));
    }

    uhd::time_spec_t time_spec = uhd::time_spec_t(usrp->get_time_now()+uhd::time_spec_t(settling));

    // start transmit worker thread
    void(*send_func)(uhd::usrp::multi_usrp::sptr, uhd::tx_streamer::sptr, const std::string&, size_t, uhd::time_spec_t);

    if (type == "double")
        send_func = &send_from_file<std::complex<double>>;
    else if (type == "float")
        send_func = &send_from_file<std::complex<float>>;
    else if (type == "short")
        send_func = &send_from_file<std::complex<short>>;
    else
        throw std::runtime_error("Unknown type " + type);

    std::thread transmit_thread([&]() {
        send_func(usrp, tx_stream, tx_file, spb, time_spec);
    });

    // recv to file
    if (type == "double")
        recv_to_file<std::complex<double>>(
            usrp, "fc64", otw, rx_file, spb, total_num_samps, rx_channel_nums, settling, time_spec);
    else if (type == "float")
        recv_to_file<std::complex<float>>(
            usrp, "fc32", otw, rx_file, spb, total_num_samps, rx_channel_nums, settling, time_spec);
    else if (type == "short")
        recv_to_file<std::complex<short>>(
            usrp, "sc16", otw, rx_file, spb, total_num_samps, rx_channel_nums, settling, time_spec);
    else {
        // clean up transmit worker
        stop_signal_called = true;
        transmit_thread.join();
        throw std::runtime_error("Unknown type " + type);
    }

    // clean up transmit worker
    stop_signal_called = true;
    transmit_thread.join();

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
