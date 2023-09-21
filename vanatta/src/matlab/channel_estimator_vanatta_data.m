fs = 2e5;
fc = 18.5e3;
fb = 500;
c = 1500;
wc = 2*pi*fc;

init_delay = 20e-3;%30e-3;
sample_offset = 24;

dd = 3;
du = 2;
d0 = dd-du;

% A0 is direct path mag, Ad is downlink path mag, Au is uplink path mag
A0 = 0.9;
Ad = 0.1;
Au = 0.2;

tau_0 = d0 / c;
tau_d = dd / c;
tau_u = du / c;

fm0_samp = fs/fb;

preamble = [0 0 1 1 1 0 1 0];
%preamble = [0 0 1 1 1 0 1 0 0 1];
%preamble = [0 0 1 1 1 0 1 0 0 1 0 0 0 1 1 1];
expected_preamble = generate_fm0_sig2(preamble,fm0_samp);
    
N_preamble_bits = length(preamble);
N_data_bits = 16;

% number of packets to decode
N_packets = 625;
packet_delay = 0; % delay in between each packet
N_tot_data_bits = N_data_bits*N_packets;
N_tot_bits = (N_data_bits+N_preamble_bits)*N_packets;

expected_data = real(read_complex_binary(strrep('../../../tx_outputs/data_prbs_order=15_len=?_packets=625.dat','?',num2str(N_data_bits))))';

%%%%% RESHAPE EXPECTED DATA INTO FORMAT FOR DFE %%%%%%
expected_data_packets = zeros(1,N_tot_bits);
% 1 2 3 4 5 6 7 8 17 18 19 ... 
preamble_index = mod(0:(N_tot_bits-1),N_preamble_bits+N_data_bits)<N_preamble_bits;
expected_data_packets(logical(preamble_index)) = repmat(preamble,1,N_packets);

data_index = (mod(0:(N_tot_bits-1),N_preamble_bits+N_data_bits)>=N_preamble_bits) ...
            .* (mod(0:(N_tot_bits-1),N_preamble_bits+N_data_bits)<(N_data_bits+N_preamble_bits));
expected_data_packets(logical(data_index)) = expected_data;
dfe_expected_data = reshape(expected_data_packets,N_preamble_bits+N_data_bits,N_packets)';
%%%%% END DFE DATA RESHAPE %%%%

% create full expected data signal
expected_data_signal = generate_fm0_sig2(expected_data_packets,fm0_samp);
expected_data_signal = expected_data_signal(1:(N_data_bits+N_preamble_bits)*fm0_samp*N_packets);

%expected_data = repmat(preamble,1,4*N_packets);

% highpass filter cutoffs
fsb1 = fb/100;
fpb1 = fb/2;
dec_fac = 2; % decimation factor before lowpass and downconversion
dfac = 5;   % donwsampling factor

% lowpass filter cutoffs
fpb1_lp = 5*fb;
fsb1_lp = 7*fb;

% % highpass for after downsampling
hpFilt = designfilt('highpassfir','PassbandFrequency',fpb1*2/(fs/(dfac*dec_fac)) ...
                    ,'StopbandFrequency',fsb1*2/(fs/(dfac*dec_fac)),'StopbandAttenuation',200,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

% lowpass after downconversion, before downsampling
lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',fpb1_lp*2/(fs/dec_fac)...
                    ,'StopbandFrequency',fsb1_lp*2/(fs/dec_fac),'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

[gdlp,w] = grpdelay(lpFilt);
gdlp = mean(gdlp);

[gdhp,w] = grpdelay(hpFilt);
gdhp = mean(gdhp);

angles = -90;
Nang = length(angles);

%%%% END DESIGN PARAMETERS %%%%

%%%% PROGRAM OPTIONS %%%%
VERBOSE = 1;
DO_PLOTS = 1;
USE_PLL = 0;
TX_LO = 0;
PURUI_PLATFORM = 0;
%%%% END PROGRAM OPTIONS %%%%

% JACK'S H, SNR, NOISE
h_median_arr = zeros(Nang,1);
h_median_snr_arr = zeros(Nang,1);
noise_median_arr = zeros(Nang,1);


% DFE H, SNR, NOISE (PRE/POST)
snr_median_pre_dfe_arr = zeros(Nang,1);
snr_median_post_dfe_arr = zeros(Nang,1);

h_median_pre_dfe_arr = zeros(Nang,1);
h_median_post_dfe_arr = zeros(Nang,1);

noise_median_pre_dfe_arr = zeros(Nang,1);
noise_median_post_dfe_arr = zeros(Nang,1);


BER = zeros(Nang,1);
BER_DFE = zeros(Nang,1);

root = '../../../rx_outputs/River PAB2 Proposal Experiments 11-09-2022/';

expected_preamble = filtfilt(lpFilt,expected_preamble')';
expected_preamble = downsample(expected_preamble,dfac*dec_fac);
expected_data_signal = downsample(expected_data_signal,dfac*dec_fac);

for n=1:Nang
    tic

    ang = angles(n);
    
    if VERBOSE
        disp("angle=");
        disp(ang);
    end

    ang_str = num2str(ang);
%     if ang >= 0
%         ang_str = strcat("+",ang_str);
%     end

    if rem(ang,1) ~= 0
        ang_str = strrep(ang_str,".",",");
    end

    if rem(ang,1) == 0
        ang_str = strcat(ang_str,',0');
    end
    
    
    filename = 'rx_single_chest_diff_004B_alone_pab2_txfmr_?deg_nicktb_18,5kfc_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_2,5m_depth_010A_purui_tx_2m_1m_hphydro_61Vrms_0.dat';

    %filename = 'rx_single_chest_pab_010B_7cm_sp_ind1,5m_+0deg_mosfet_18,5kfc_siggen_data_1kbps_usrp_2,5m_depth_3m_u2b_0,5m_hphydro_0.dat';
    filepath = strcat(root,strrep(filename,'?',ang_str));
    
    if TX_LO
        filepath = strrep(filepath,'.dat','.00.dat');
    end

    yr = read_complex_binary(filepath);
    
    if PURUI_PLATFORM
        rx_signals = real(yr).';
    else
        sig = yr(sample_offset:end).';
        rx_signals = real(sig)-imag(sig);
    end

    rx_signals = decimate(rx_signals,dec_fac);
    fs_n = fs/dec_fac;
    rx_len = length(rx_signals);
    
    %%%% CARRIER FREQUENCY AND PHASE EXTRACTION %%%%
    % have had some issues with it in the past and since RX and TX USRPs are 
    % synchronized in most experiments directly using the known fc works fine
    
    Nfft = 100*fs_n;
    rx_fft = fft(rx_signals',Nfft)';
    fft_mag = abs(rx_fft);
    max_search = [round(Nfft/fs_n*(fc-1)):round(Nfft/fs_n*(fc+1))];
    [maxval,mindex] = max(fft_mag(:,max_search),[],2); % max in each row
    carrier_phase = angle(rx_fft(max_search(mindex)'));
    carrier_freq = fs_n/Nfft*max_search(mindex)';
    
    %carrier_freq = fc;
    %carrier_phase = 0;
    
    % generate the time series and local oscillator
    t = [0:1/fs_n:(rx_len-1)/fs_n];
    lo = exp(1j*(2*pi*carrier_freq*t+carrier_phase));
    
    if TX_LO
        lo = real(read_complex_binary(strrep(filepath,'.00.dat','.01.dat')))';
        lo = lo(sample_offset:end);
          

        lo = decimate(lo,dec_fac);

%         corr = xcorr(rx_signals',lo');
%         corr = corr(length(rx_signals):end);
%         [max_val,lo_start] = max(corr);
% 
%         rx_signals = rx_signals(lo_start:end);
%         lo = lo(1:end-lo_start+1);
        
    %%%% SOFTWARE PLL %%%%
    % this block replaces t and lo if USE_PLL = 1
    elseif USE_PLL
        t_tot = rx_len/fs;
        
        ph = zeros(1,rx_len);
        ph_est = zeros(1,rx_len);
        lp = zeros(1,rx_len);
        y = zeros(1,rx_len);
        integ = zeros(1,rx_len);
        
        Bn = 1e-2*fs;
        damp = 1/sqrt(2);
        
        k0 = 1;
        kd = 0.5;
        kp = 1/(kd*k0)*4*damp/(damp+1/(4*damp))*Bn/fs;
        ki = 1/(kd*k0)*4/(damp+1/(4*damp))^2*(Bn/fs)^2;
        
        integ_out = 0;
        ph_est(1) = carrier_phase;
        
        for i = 1:rx_len-1
            t(i) = t_tot*i/rx_len;
            % input signal
            y(i) = rx_signals(1,i);%sin(2*pi*fc*t(n)+pi);
        
            % phase detect
            ph(i) = kd*y(i)*imag(lo(i));
        
            % loop filter
            integ_out = ki*ph(i)+integ_out;
            lp(i) = kp*ph(i) + integ_out;
        
            % vco
            ph_est(i+1) = ph_est(i) + k0*lp(i);
            lo(i+1) = exp(-1j*(2*pi*carrier_freq*t_tot*(i+1)/rx_len+ph_est(i)));
        end
    
        t(end) = t_tot;
        
        figure(10);
        plot(t,rx_signals(1,:));
        hold on;
        plot(t,real(lo));
    end

    % downconvert
    rx_baseband = rx_signals.*lo;
    
    % lowpass filtering both removes the 2fc term and anti-alias filters
    % filtfilt used for 0 group delay filtering
    rx_baseband = fftfilt(lpFilt,rx_baseband')';
    rx_baseband = downsample(rx_baseband,dfac);

    fs_n = fs_n/dfac;
    fm0_samp = fs_n/fb;
    
    rx_baseband = fftfilt(hpFilt,rx_baseband')';
    rx_baseband = rx_baseband(round(init_delay*fs_n+gdlp+gdhp-fm0_samp):end);
    
    sig_sec = rx_baseband;
%     Nfft = 2^nextpow2(length(sig_sec));
%     fft_sig = fft(sig_sec,Nfft);
%     [subcar_peak,subcar_peak_index] = max(fft_sig);
%     subcar_peak_f = subcar_peak_index/Nfft*fs_n;
%     subcar_peak_phase = angle(subcar_peak);
%     f = [-fs_n/2:fs_n/Nfft:fs_n/2-fs_n/Nfft];
    
    t_window = 0.1;
    window = chebwin(t_window*fs_n);
    Nfft = 2^nextpow2(length(window));
    
    if DO_PLOTS
        figure(1);
        hold on;
        [pxx,f] = pwelch(sig_sec,window,[],Nfft,fs_n);
        plot(f,10*log10(pxx));
        grid on;
        grid minor;

    
        % disp("Total Average Power (dB):");
        % disp(num2str(10*log10(tot_pow)));
        
        figure(2);
        subplot(2,1,1);
        hold on;
        title("Real Part RX Baseband")
        plot(real(sig_sec));
        ylabel("Amplitude (V)");
        subplot(2,1,2);
        hold on;
        plot(imag(sig_sec));
        ylabel("Amplitude (V)");
        xlabel("Samples");
        title("Imag Part RX Baseband");
    end
    
    
%% CORRELATE, ESTIMATE, DECODE, and COMPUTE BER %%
    
    % N_data_bits = 0;
    % 
    % % resample for even integer fm0_samp value
    % % with fs_n/fb divisble by 2 this code is not run
    % if mod(fs_n/fb,2) ~= 0
    %     p = (ceil(fs_n/fb)+1)*fb;
    %     q = fs_n;
    %     
    %     rx_baseband = resample(rx_baseband,p,q);
    %     
    %     fs_n = p;
    % end
    
    % % estimates after channel projection for individual packets
    preamble_len = fm0_samp*N_preamble_bits;
    data_len = fm0_samp*N_data_bits;                 % known data length in samples
    packet_len = data_len+preamble_len;     % packet length in samples
   
    data_estimates = zeros(1,data_len*N_packets);
    decoded_data = zeros(1,N_packets*N_data_bits);
    packet_estimates = zeros(1,packet_len*N_packets);
    
    channel_estimates = zeros(1,N_packets);
    channel_snrs = zeros(1,N_packets);
    noise_power = zeros(1,N_packets);
    
    decode_preamble = expected_preamble-mean(expected_preamble);
    
    packet_delay_adj = floor(packet_delay*fs_n);

    % tx norm is length of preamble for binary keying (-1,+1)
    tx_norm = sum(abs(decode_preamble).^2);

    fm0_half_samp = ceil(fm0_samp/2);

    %%% full data sequence correlation to find global preamble start %%%
    rcorr = xcorr(real(rx_baseband)',expected_data_signal');
    icorr = xcorr(imag(rx_baseband)',expected_data_signal');

    corr_tot = rcorr(length(rx_baseband):end)+1j*icorr(length(rx_baseband):end);
    abs_corr = abs(corr_tot);
    [preamble_max,global_preamble_start] = max(abs_corr);
    abs_corr = abs_corr/preamble_max;
    %%%% end full data sequence correlation %%%%
    
    %global_preamble_start = global_preamble_start + fm0_samp/2;
    

    if DO_PLOTS
        figure(3);
        plot(abs_corr/100);
        hold on;
        plot(real(rx_baseband));
        plot([zeros(1,global_preamble_start) expected_data_signal/1000]);
    end
    
    
    % cutoff rx_baseband where it begins
    rx_baseband = rx_baseband(global_preamble_start-fm0_samp:end);
        
    if DO_PLOTS
        figure(4);
    end
    % CORRELATION AND DECODING %
    for pnum=1:N_packets
        % remove +sample_delay_adj*(pnum-1) for single correlation
        begdex = (pnum-1)*packet_len+1+packet_delay_adj*(pnum-1);
        endex = pnum*packet_len+packet_delay_adj*(pnum-1);
        end_preamble_dex = (pnum-1)*packet_len+floor(1.45*preamble_len)+packet_delay_adj*(pnum-1);

        beg_data_dex = (pnum-1)*data_len+1;
        end_data_dex = pnum*data_len+fm0_samp;
    
        beg_bit_dex = (pnum-1)*N_data_bits+1;
        end_bit_dex = pnum*N_data_bits;
        
        % perform cross correlation for packet start
%         [rcorr,rlags] = xcorr(real(rx_baseband(begdex:end_preamble_dex)).',decode_preamble.');
%         [icorr,ilags] = xcorr(imag(rx_baseband(begdex:end_preamble_dex)).',decode_preamble.');
%         % removes tails of correlation 
%         corr_tot = rcorr(end_preamble_dex-begdex+1:end)+1j*icorr(end_preamble_dex-begdex+1:end);
%         abs_corr = abs(corr_tot);
%         % find maximum correlation and begin decoding from there
%         [preamble_max,preamble_start] = max(abs_corr); 
        
        %%% comment out for correlation at every packet
        preamble_start = fm0_samp;

        if DO_PLOTS
            clf;
            
            plot(abs_corr/30);
            hold on;
            plot(10*real(rx_baseband(begdex:end_preamble_dex)));
            plot([zeros(1,preamble_start) decode_preamble/100]);
        end
        %xlim([0 1000]);

%         if pnum == 207
%             test = 1;
%         end

        % slice out data packet found from correlation
        packet = rx_baseband(begdex+preamble_start-fm0_half_samp:endex+preamble_start+fm0_half_samp);
        % remove the mean
        packet = packet - mean(packet(1:preamble_len));
        % slice out preamble
        packet_preamble = packet(fm0_half_samp+1:fm0_half_samp+preamble_len);
        % slice out data
        packet_data = packet(preamble_len+1:end);
        
%         %remove edges
%         packet_preamble = remove_edges(packet_preamble,preamble1,0.05,fb,fs_n);
%         decode_preamble = remove_edges(decode_preamble,preamble1,0.05,fb,fs_n);
%         % recompute tx_norm since some samples have been removed
%         tx_norm = sum(abs(decode_preamble).^2);

        % compute the channel estimate 
        channel_estimates(pnum) = sum(packet_preamble.*conj(decode_preamble))/tx_norm;

        % SNR calculation 
        noise_est = packet_preamble-channel_estimates(pnum).*decode_preamble;
        noise_est = reshape(noise_est,[fm0_samp N_preamble_bits]);
        
        noise_est_per_bit = mean(noise_est);
        noise_power(pnum) = var(noise_est_per_bit);

        channel_snrs(pnum) = abs(channel_estimates(pnum))^2/noise_power(pnum);
        
        % extract estimate of data only
        data_estimates(beg_data_dex:end_data_dex) = packet_data.*conj(channel_estimates(pnum))./abs(channel_estimates(pnum)).^2;
        packet_estimates(begdex:endex) = packet(fm0_half_samp+1:end-fm0_half_samp).*conj(channel_estimates(pnum))./abs(channel_estimates(pnum)).^2;   

        % decode packet
        if (mod(fm0_samp,2) == 0)
            % should always enter this branch, camera_decode written by
            % saad and waleed
            bits = camera_decode(data_estimates(beg_data_dex:end_data_dex).',fm0_samp,N_data_bits);
        else
            error("fm0_samp is not divisible by 2");
        end
        
        % place bits into decoded data matrix
        decoded_data(beg_bit_dex:end_bit_dex) = bits;
    end
    
    if DO_PLOTS
        figure(5);
        hold on;
        plot(channel_estimates,'x','linewidth',2);
        
        xlim([-1 1]*1e-3);
        ylim([-1 1]*1e-3);
        
        xlabel("Real");
        ylabel("Imaginary");
        
        %axis equal;
        grid on;
        grid minor;
        title("Channel Estimates Constellation");
    end

    h_median_arr(n) = median(channel_estimates);
    h_median_db = 20*log10(h_median_arr(n));
    h_median_snr_arr(n) = 10*log10(median(channel_snrs));
    noise_median_arr(n) = 10*log10(median(noise_power));
    
    min_BER = 1/N_tot_data_bits;
    BER(n) = sum(decoded_data ~= expected_data)/(N_data_bits*N_packets);
    BER(BER == 0) = min_BER;

    rx_baseband = rx_baseband(fm0_samp+1:N_tot_bits*fm0_samp+fm0_samp);
    %%
    N_training = 300;
    [weights,ber_fin_b,ber_fin_a,snr_pre_final,snr_post_final,h_mag_pre_final,h_mag_post_final,noise_pre_final,noise_post_final] = ... 
        DFE_500_vanatta(packet_estimates.',dfe_expected_data,N_training,N_packets-N_training,1,0,1,fb,fs_n,1);
    
%     N_training = 1;
%     [weights,ber_fin_b,ber_fin_a,snr_pre_final,snr_post_final,h_mag_pre_final,h_mag_post_final,noise_pre_final,noise_post_final] = ...
%         DFE_500_vanatta(packet_estimates.',dfe_expected_data,N_training,N_packets-N_training,1,weights,1,fb,fs_n,1);
    
    BER_DFE(n) = ber_fin_a;

    snr_median_pre_dfe_arr(n) = 10*log10(snr_pre_final);
    snr_median_post_dfe_arr(n) = 10*log10(snr_post_final);
    
    h_median_pre_dfe_arr(n) = 10*log10(h_mag_pre_final);
    h_median_post_dfe_arr(n) = 10*log10(h_mag_post_final);

    noise_median_pre_dfe_arr(n) = 10*log10(noise_pre_final);
    noise_median_post_dfe_arr(n) = 10*log10(noise_post_final);

    toc
end
%% PLOT VS ANGLE
if length(angles) > 1
%     figure(6);
%     hold on;
%     plot(angles,20*log10(abs(h_median_arr)));
%     grid on;
%     grid minor;
%     xlabel("Angle (deg)");
%     ylabel("Channel Mag (dB20)");
%     
%     figure(7);
%     hold on;
%     plot(angles,h_median_snr_arr);
% %     plot(angles,20*log10(abs(h_median_arr))-noise_median_arr);
%     grid on;
%     grid minor;
%     xlabel("Angle (deg)");
%     ylabel("Channel SNR (dB20)");
% 
%     figure(8);
%     semilogy(h_median_snr_arr,BER,'o','LineWidth',2);
%     hold on;
%     grid on;
%     grid minor;
%     xlabel("Median Bit SNR (dB)");
%     ylabel("Total BER");
% 
%     figure(9);
%     hold on;
%     plot(angles,noise_median_arr);
%     grid on;
%     grid minor;
%     title("Median Noise Power (dB) vs. Angle (deg)");
%     xlabel("Angle (deg)");
%     ylabel("Median Noise Power (dB)");

    figure(6);
    
    plot(angles,h_median_pre_dfe_arr);
    grid on;
    grid minor;
    title("Median Channel Magnitude PRE DFE");
    xlabel("Angle (deg)");
    ylabel("Magnitude (dB)");
    hold on;
    
    figure(7);
    plot(angles,snr_median_pre_dfe_arr);
%     plot(angles,20*log10(abs(h_median_arr))-noise_median_arr);
    grid on;
    grid minor;
    xlabel("Angle (deg)");
    ylabel("Magnitude (dB)");
    title("Median SNR PRE DFE");
    hold on;
% 
%     figure(8);
%     semilogy(h_median_snr_arr,BER,'o','LineWidth',2);
%     hold on;
%     grid on;
%     grid minor;
%     xlabel("Median Bit SNR (dB)");
%     ylabel("Total BER");

    figure(9);
    plot(angles,noise_median_pre_dfe_arr);
    grid on;
    grid minor;
    title("Median Noise Power PRE DFE");
    xlabel("Angle (deg)");
    ylabel("Magnitude (dB)");
    hold on;
end
% (abs(comb_ch_est(2:end))-abs(n1_ch_est+n2_ch_est(2:end)))./abs(comb_ch_est(2:end))


%% EXPORT DATA %%
% t = [0:1/fs:1-1/fs];
% t0 = t - init_delay;
% 
% len_packet1 = length(expected_preamble1)*n_data_reps1/fs;
% len_packet2 = length(expected_preamble2)*n_data_reps2/fs;
% 
% do_vanatta = 0;
% 
% data_real = zeros(1,length(t),'like',t);
% data_imag = zeros(1,length(t),'like',t);
% 
% for n=1:N_trials
%     seg_ch1 = data(t0-(n-1)*len_packet1-2*(n-1)*len_packet2-(3*n-3)*packet_delay,preamble1,fb,n_data_reps1);
%     seg_ch2 = data(t0-n*len_packet1-2*(n-1)*len_packet2-(3*n-2)*packet_delay,preamble2,fb,n_data_reps2);
%     seg_comb = data(t0-n*len_packet1-(2*n-1)*len_packet2-(3*n-1)*packet_delay,preamble2,fb,n_data_reps2);
% 
%     data_real = data_real+seg_ch1+seg_comb/(do_vanatta+1);
%     data_imag = data_imag+seg_ch2+seg_comb/(do_vanatta+1);
% end
% 
% if do_vanatta
%     data_tot = (1+1j)*(data_real+data_imag)*0.7;
%     data_tot(real(data_tot)+imag(data_tot)<0) = 0;
%     data_tot(logical((t < init_delay)+(t > init_delay+(3*N_trials-1)*packet_delay+N_trials*(len_packet1+2*len_packet2)))) = 0;
%     write_complex_binary(data_tot,strrep("../../tx_outputs/vanatta_channel_estimating_data_?bps.dat","?",num2str(round(fb))));
% else
%     data_tot = ((data_real)+1j*(data_imag))*0.7;
%     data_tot(real(data_tot)+imag(data_tot)<0) = 0;
%     data_tot(logical((t < init_delay)+(t > init_delay+(3*N_trials-1)*packet_delay+N_trials*(len_packet1+2*len_packet2)))) = 0;
%     write_complex_binary(data_tot,strrep("../../tx_outputs/array_channel_estimating_data_?bps.dat","?",num2str(round(fb))));
% end
% % 
% figure(1);
% hold on;
% plot(t,real(data_tot));
% plot(t,imag(data_tot));

% out = remove_edges(expected_preamble1,preamble1,0.05,fb,fs);
% plot(out);
% hold on;
% plot(expected_preamble1);

function out = data(t,code,fb,nreps)
    fm0_code = generate_fm0_sig2(code,2);
    code_len = length(fm0_code);
    out = zeros(size(t),'like',t);

    for n=1:nreps*code_len
        out(logical((t < n*1/(2*fb)).*(t >= (n-1)*1/(2*fb)))) = fm0_code(mod(n-1,code_len)+1);
    end
end

% removal is a percentage of the bit period
function out = remove_edges(data,expected_code,percent_removal,fb,fs)
    if length(data) ~= length(expected_code)*fs/fb
        error("Data length is not equal to expected length based on fb and fs given")
    end

%     fm0_code = generate_fm0_sig2(expected_code,2);
    
    % array of how many expected edges there are and their indexes (without
    % beginning edge)
    edge_indexes = zeros(length(expected_code)+length(expected_code(expected_code==0)),1);
    
    nedge = 1;
    for n=1:length(expected_code)
        % fill edge center indexes
        bit = expected_code(n);
        
        if bit
            edge_indexes(nedge) = n*fs/fb;
            nedge = nedge + 1;
        else
            edge_indexes(nedge) = n*fs/fb-fs/(2*fb);
            nedge = nedge + 1;
            edge_indexes(nedge) = n*fs/fb;
            nedge = nedge + 1;
        end
    end

    sample_expansion = 2*floor(fs/fb*percent_removal/2)+1;
    samples_to_remove = repelem(edge_indexes,sample_expansion);
    samples_to_remove = samples_to_remove + repmat([-(sample_expansion-1)/2:1:(sample_expansion-1)/2]',length(edge_indexes),1);
    samples_to_remove(end-(sample_expansion-1)/2+1:end) = [];
    samples_to_remove = cat(1,[1:(sample_expansion-1)/2]',samples_to_remove);

    data(samples_to_remove) = [];
    out = data;
end

function y = sinc_interp(x,s,u)
    % Interpolates x sampled sampled at "s" instants
    % Output y is sampled at "u" instants ("u" for "upsampled")
    % (EXPECTS x, s, and u to be ROW VECTORS!!)

    % Find the period of the undersampled signal
    T = s(2)-s(1);

    % When generating this matrix, remember that "s" and "u" are
    % passed as ROW vectors and "y" is expected to also be a ROW
    % vector. If everything were column vectors, we'd do.
    %
    % sincM = repmat( u, 1, length(s) ) - repmat( s', length(u), 1 );
    %
    % So that the matrix would be longer than it is wide.
    % Here, we generate the transpose of that matrix.
    sincM = repmat( u, length(s), 1 ) - repmat( s', 1, length(u) );

    % Equivalent to column vector math:
    % y = sinc( sincM'/T )*x';
    y = x*sinc( sincM/T );
end
