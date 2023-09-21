%%%% DESIGN PARAMETERS %%%%
Nel = 1;
fs = 2e5;               % USRP sampling rate (is changed throughout script due to downsampling)
fc = 18.5e3;              % carrier frequency
fb = 1e3;

% generate low and highpass filters
% lowpass is performed before downsampling as an anti-aliasing filter
% highpass is performed after downsampling as it is a higher order filter
fsb1 = fb/10;
fpb1 = fb;
dfac = 5;   % donwsampling factor

% highpass for after downsampling
hpFilt = designfilt('highpassfir','PassbandFrequency',fpb1*2/(fs/dfac) ...
                    ,'StopbandFrequency',fsb1*2/(fs/dfac),'StopbandAttenuation',40,'PassbandRipple',0.1);

% lowpass after downconversion, before downsampling
lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',fs/dfac/8*2/fs...
                    ,'StopbandFrequency',fs/dfac/8*3/2*2/fs,'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1);

%%%% END DESIGN PARAMETERS %%%%

% the root of the filename of the rx data. Remove the ending _0 _1 _2 from
% filename and place here
folder = '~/Documents/MIT/sk/oceans/vanatta/rx_outputs/River Switch Van Atta Tests 06-20-2022/';
file = 'rx_backscatter_vanatta_switch_TS5A_pab_003A_006A_ind_0deg_18,5kfc_1kmod_2m_depth_0.dat';
root = strcat(folder,file);

% initializes size
sig = read_complex_binary(strcat(root));
%sig_synth = read_complex_binary(root_synth);

sig = sig(24:end);

rx_len = length(sig);
% Nel x rx_len size matrix of input signals, where each row is time-series on an individual array element
rx_signals = zeros(Nel,rx_len); 

% synthetic data
% rx_signals = zeros(Nel,length(sig_synth));
% rx_signals(1,:) = awgn(sig_synth,20,'measured');

% reads files of root + _0.dat, root + _1.dat, and root + _2.dat, depending
% on Nel
if Nel > 1
    for fnum=0:Nel/2-1
        sig = read_complex_binary(strcat(root,'_',int2str(fnum),'.dat'));
        
        rx_signals(2*fnum+1,:) = real(sig(1:rx_len));
        rx_signals(2*fnum+2,:) = imag(sig(1:rx_len));
    end
else
%     sig = read_complex_binary(strcat(root));
    rx_signals(1,:) = real(sig);
end

rx_len = length(rx_signals(1,:));

%%%% CARRIER FREQUENCY AND PHASE EXTRACTION %%%%
% have had some issues with it in the past and since RX and TX USRPs are 
% synchronized in most experiments directly using the known fc works fine

%Nfft = 100*fs;
% rx_fft = fft(rx_signals',Nfft)';
% fft_mag = abs(rx_fft);
% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(fft_mag(:,max_search),[],2); % max in each row
% carrier_phase = angle(rx_fft(max_search(mindex)'));
% carrier_freq = fs/Nfft*max_search(mindex)';

clear rx_fft fft_mag;
carrier_freq = fc;
carrier_phase = 0;

%%%% OVERLAP ADD METHOD FOR FILTERING LONG INPUT SIGNALS %%%%
% if you want to learn more about this method there is a section on it in
% Discrete Time Signal Processing by Oppenheim & Schafer (Jack has book)

L = 1e6;    % window length
P = length(lpFilt.Coefficients);    % filter order (filter length)

rx_baseband = zeros(Nel,rx_len+P-1);

% inside for loop rx_signals is multiplied by a local oscillator and then
% fft-filtered by lowpass generated in design parameters section 
for seg = 1:ceil(length(rx_signals(1,:))/L)
    endex = seg*L;
    if endex > length(rx_signals(1,:))
        endex = length(rx_signals(1,:));
    end

    rx_segment = padarray(rx_signals(:,(seg-1)*L+1:endex),[0 P-1],0,'post');

    seg_len = length(rx_segment(1,:));
    
    % generate the time series and local oscillator
    t = [(seg-1)*L/fs:1/fs:(seg_len+(seg-1)*L-1)/fs];
    lo = exp(1j*(2*pi*carrier_freq*t+carrier_phase));
    % factor of 2 comes from cosine expansion
    rx_baseband_seg = 2*rx_segment.*lo;
    
    % lowpass filtering both removes the 2fc term and anti-alias filters
    % the signal to prepare for downsampling
    rx_baseband_seg = filtfilt(lpFilt,rx_baseband_seg')';
    
    %expected_preamble = repelem(preamble,fm0_samp/2);
    
    if seg == 1
        rx_baseband(:,1:L+P-1) = rx_baseband_seg;
        %rx_baseband(:,(seg-1)*L+1:seg*L+P-2+1) = rx_baseband+
    else
        % overlap and add
        overlap_begdex = (seg-1)*L+1;
        overlap_endex = (seg-1)*L+P-1;

        remain_begdex = (seg-1)*L+P;
        remain_endex = seg*L+P-1;

        if remain_endex > rx_len+P-1
            remain_endex = rx_len+P-1;
        end

        rx_baseband(:,overlap_begdex:overlap_endex) = rx_baseband(:,overlap_begdex:overlap_endex)+rx_baseband_seg(:,1:P-1);
        rx_baseband(:,remain_begdex:remain_endex) = rx_baseband_seg(:,P:end);
    end
end

% expected_preamble = filtfilt(lpFilt, expected_preamble);

% downsample by 5
rx_baseband=downsample(rx_baseband',dfac)';
% expected_preamble = downsample(expected_preamble,dfac);
fs = fs/dfac;   % change sampling rate to reflect downsampled data

% highpass filter at lower sampling rate
rx_baseband = filtfilt(hpFilt,rx_baseband')';
%expected_preamble = fftfilt(hpFilt,expected_preamble);

t_window = 1;
window = chebwin(t_window*fs);
Nfft = length(window);

[pxx,f] = pwelch(rx_baseband,window,[],Nfft,fs,'power');
plot(f,10*log10(pxx));

figure;
plot(real(rx_baseband));

tot_pow = bandpower(rx_baseband);

disp("Total Average Power (dB):");
disp(num2str(10*log10(tot_pow)));

% hold on;
% plot(expected_preamble)
% plot(real(rx_baseband(1,:)));

% this is the hacky code I wrote to detect where the signal starts
% index_from and value_larger_than might need to be changed depending on
% the character of the signal 

% index_from = 1000;
% value_larger_than = 0.07;
% start_index = find(abs(rx_baseband(1,index_from:end)) > value_larger_than, 1) + index_from - 500;
% if(start_index < 1)
%     start_index = 1;
% end
% %start_index = 1;
% % start processing from "start" of data
% rx_baseband = rx_baseband(:,start_index:end);