sig_0 = read_complex_binary('../../rx_outputs/River PAB2 Van Atta 8 11-09-2022/noise_test_diff_11m_10m_usrp_18,5k_XLI_amp_61Vrms.dat');

sty = 'real';

if isequal(sty,'diff')
    sig = real(sig_0(24:end))'-imag(sig_0(24:end))';
elseif isequal(sty,'real')
    sig = real(sig_0(24:end))';
else
    sig = imag(sig_0(24:end))';
end

fs = 2e5;
dfac = 4;   % downsampling factor

sig = decimate(sig,dfac);
fs = fs/dfac;
t = [0:1/fs:length(sig)/fs-1/fs];
% figure(1);
% plot(sig);
fc = 18.5e3;

window_size = floor(length(sig)/20);
window = chebwin(window_size);
Nfft = fs*10;

[pxx,f] = pwelch(sig,window,[],Nfft,fs,'power');

figure;
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");

%% downconvert
fb = 500;

% highpass filter cutoffs
fsb1 = fb/100;
fpb1 = fb/2;

% lowpass filter cutoffs
fpb1_lp = 3*fb;
fsb1_lp = 5*fb;

% % highpass for after downsampling
hpFilt = designfilt('highpassfir','PassbandFrequency',fpb1*2/fs ...
                    ,'StopbandFrequency',fsb1*2/fs,'StopbandAttenuation',80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

% lowpass after downconversion, before downsampling
lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',fpb1_lp*2/fs...
                    ,'StopbandFrequency',fsb1_lp*2/fs,'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

init_delay = 50e-3;

[gdlp,w] = grpdelay(lpFilt);
gdlp = mean(gdlp);

[gdhp,w] = grpdelay(hpFilt);
gdhp = mean(gdhp);

%%%% CARRIER FREQUENCY AND PHASE EXTRACTION %%%%
% have had some issues with it in the past and since RX and TX USRPs are 
% synchronized in most experiments directly using the known fc works fine

Nfft = 500*fs;
rx_fft = fft(sig',Nfft)';
fft_mag = abs(rx_fft);
[maxval,mindex] = max(fft_mag); % max in each row
carrier_freq = fs/Nfft*mindex;

 % generate the time series and local oscillator
lo = exp(1j*(2*pi*carrier_freq*t));
% downconvert + filter
rx_baseband = sig.*lo;
rx_baseband = fftfilt(lpFilt,rx_baseband')';
rx_baseband = fftfilt(hpFilt,rx_baseband')';

rx_baseband = rx_baseband(floor(init_delay*fs+gdlp+gdhp-40):end);

figure;
subplot(2,1,1);
hold on;
plot(real(rx_baseband));
subplot(2,1,2);
hold on;
plot(imag(rx_baseband));

disp("Baseband RMS Noise = ");
disp(rms(rx_baseband));

disp("RMS PB = ");
disp(rms(sig));
