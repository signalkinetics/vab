
sig = real(read_complex_binary('../../rx_outputs/WHOI Experiments 11-30-2022/barge_noise_test_switching_off_amp_off_2m_2m_1m_sep_hphydro_0.dat'));
fs = 2e5;
sig = sig(24:end);

%sig = sig(16e4:end);

t = [0:1/fs:length(sig)/fs-1/fs];

window_size = floor(length(sig)/10);
window = chebwin(window_size,120);
Nfft = 2^nextpow2(window_size);

rbw = enbw(window,fs);

[pxx,f] = pwelch(sig,window,[],Nfft,fs,'power');
pxx = pxx/50;

% figure(1);
% plot(sig);

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");

% compute phase noise
offsets = [ceil(3*rbw/10)*10:20:5e3];

[carrier_val,carrier_index] = max(pxx);
carrier_freq = f(carrier_index);

phasenoise = 10*log10(pxx(round(carrier_index+offsets/fs*Nfft))/carrier_val);
% compute FM0 spectrum
fb = 500;
Nbits = 1024;

code = randi([0 1],1,Nbits); 
fm0_seq = generate_fm0_sig2(code,fs/fb);
% t = [0:1/fs:N_bits/fb-1/fs];

window_size = floor(length(fm0_seq)/10);
window = chebwin(window_size);
Nfft = 2^nextpow2(window_size);

[pxx_fm0,f] = pwelch(fm0_seq,window,[],[],fs);

figure(3);
hold on;
plot(offsets,phasenoise);
%plot(f(1:floor(offsets(end)*Nfft/fs)),10*log10(pxx_fm0(1:floor(offsets(end)*Nfft/fs))));
grid on;
grid minor;
xlabel("Freq Offset from Carrier (Hz)");
ylabel("Phase noise (dBc/Hz)");
title(strcat("Phase Noise w/ Carrier Freq = ",num2str(carrier_freq/1e3)," kHz"));




