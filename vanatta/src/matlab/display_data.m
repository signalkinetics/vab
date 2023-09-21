sig_0 = read_complex_binary('../../rx_outputs/River PAB2 Van Atta 12-20-2022/fixed_vanatta4x2_nostag_006B_006F_006A_006C_x_001A_004A_004B_004D_chest_txfmr_nicktb_siggen_18,5kfc_0,0deg_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_2,5m_depth_005B_purui_tx_60Vrms_18m_18m_1m_foam_sep_hphydro_diff_0.00.dat');
sig = real(sig_0(24:end))-imag(sig_0(24:end));
fs = 2e5;
fc = 18.5e3;

%sig = rx_baseband;
t = [0:1/fs:length(sig)/fs-1/fs];
%sig = awgn(cos(2*pi*fc*t),30);
figure(1);
hold on;
plot(sig);


window_size = floor(length(sig)/100);
window = chebwin(window_size);
Nfft = 2^nextpow2(fs*10);

[pxx,f] = pwelch(sig,window,[],Nfft,fs,'power');
% 
% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';


figure(2);

plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");
ylabel("dB re 1uPa^2/Hz");
% xlim([1e3 100e3]);
hold on;
%%
% lpFilt = designfilt('lowpassfir' ...
%                     ,'PassbandFrequency',20e3*2/(fs)...
%                     ,'StopbandFrequency',30e3*2/(fs),'StopbandAttenuation' ...
%                     ,200,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

dfac = 4;
%new_sig = fftfilt(lpFilt,sig);
new_sig=sig;
new_sig = decimate(new_sig,dfac);
fs = fs/dfac;

window_size = floor(length(new_sig)/1);
window = chebwin(window_size);
Nfft = 2^nextpow2(fs*10);

[pxx,f] = pwelch(new_sig,window,[],Nfft,fs);

% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");


