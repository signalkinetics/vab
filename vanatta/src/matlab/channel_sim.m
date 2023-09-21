fs = 2e5;
fc = 18.5e3;
fb = 1e3;
c = 1500;
wc = 2*pi*fc;

%theta = 0; % array rotation
if theta > 90
    return
end

ds = 0.07; % element spacing
delta_d = ds*sin(theta/180*pi);

dd = 3;
du = 2;

dd1 = dd+delta_d/2;
dd2 = dd-delta_d/2;

du1 = du+delta_d/2;
du2 = du-delta_d/2;

d0 = 1;

% A0 is direct path mag, Ad is downlink path mag, Au is uplink path mag
A0 = 0.9;
Ad = 0.1;
Au = 0.15;

tau_0 = d0 / c;
tau_d1 = dd1 / c;
tau_d2 = dd2 / c;

tau_u1 = du1 / c;
tau_u2 = du2 / c;

Nel = 1;

preamble1 = [0 0 1 1 1 0 1 0];
expected_preamble1 = generate_fm0_sig2(preamble1,fs/fb);

preamble2 = preamble1;
expected_preamble2 = generate_fm0_sig2(preamble2,fs/fb);

expected_preamble_comb = expected_preamble1;

% m = 0.1;
t_len = 1;
t = [0:1/fs:t_len-1/fs];

%plot(t,data(t-0.1,fb,1));
n_data_reps1 = 5;
n_data_reps2 = 5;
data1 = data(t-tau_u1,preamble1,fb,n_data_reps1);
data2 = data(t-(length(expected_preamble1)*n_data_reps1/fs+1e-3)-tau_u2,preamble2,fb,n_data_reps2);

data_comb_part1 = data(t-(length(expected_preamble1)*n_data_reps1/fs+length(expected_preamble2)*n_data_reps2/fs+2e-3)-tau_u1,...
                    preamble1,fb,n_data_reps1);
data_comb_part2 = data(t-(length(expected_preamble1)*n_data_reps1/fs+length(expected_preamble2)*n_data_reps2/fs+2e-3)-tau_u2,...
                    preamble2,fb,n_data_reps2);

%data_comb = data_comb_part1+data_comb_part2;

yr = A0*cos(wc*(t-tau_0)).*heaviside(t-tau_0) + ...
        Ad*Au*data1.*cos(wc*(t-tau_d1-tau_u1)).*heaviside(t-tau_d1-tau_u1) + ...
        Ad*Au*data2.*cos(wc*(t-tau_d2-tau_u2)).*heaviside(t-tau_d2-tau_u2) + ...
        Ad*Au*data_comb_part1.*cos(wc*(t-tau_d1-tau_u1)).*heaviside(t-tau_d1-tau_u1) + ...
        Ad*Au*data_comb_part2.*cos(wc*(t-tau_d2-tau_u2)).*heaviside(t-tau_d2-tau_u2);
%plot(t,yr);
% carrier = 0;
% pb_sig = (1+m*data).*carrier;
% pb_sig = pb_sig';

load 3m_pb_ch.mat;
channel = ch_pb2;
% figure;
% periodogram(channel);

% out = conv(pb_sig,channel);
% 
% t_out = [0:1/fs:length(out)/fs-1/fs];
% downconv_out = fftfilt(lpFilt, cos(2*pi*fc*t_out)'.*out);

% figure;
% subplot(2,1,1);
% plot(pb_sig);
% ylim([-2 2]);
% subplot(2,1,2);
% plot(downconv_out);

% t_window = 0.1;
% window = chebwin(t_window*fs);
% Nfft = 2^nextpow2(length(window));

% figure;
% hold on;
% [pxx,f] = pwelch(data,window,[],Nfft,fs,'power');
% plot(f,10*log10(pxx));
% grid on;
% grid minor;

% generate low and highpass filters
% lowpass is performed before downsampling as an anti-aliasing filter
% highpass is performed after downsampling as it is a higher order filter
fsb1 = fb/100;
fpb1 = fb/4;
dfac = 1;   % donwsampling factor

fpb1_lp = 18e3;
fsb1_lp = 20e3;

% % highpass for after downsampling
% hpFilt = designfilt('highpassfir','PassbandFrequency',fpb1*2/(fs/dfac) ...
%                     ,'StopbandFrequency',fsb1*2/(fs/dfac),'StopbandAttenuation',80,'PassbandRipple',1,'DesignMethod','kaiserwin');

% lowpass after downconversion, before downsampling
lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',fpb1_lp*2/fs...
                    ,'StopbandFrequency',fsb1_lp*2/fs,'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

%%%% END DESIGN PARAMETERS %%%%

sig = yr;

rx_len = length(sig);
% Nel x rx_len size matrix of input signals, where each row is time-series on an individual array element
rx_signals = zeros(Nel,rx_len);
rx_signals(1,:) = real(sig);

%%%% CARRIER FREQUENCY AND PHASE EXTRACTION %%%%
% have had some issues with it in the past and since RX and TX USRPs are 
% synchronized in most experiments directly using the known fc works fine

Nfft = 100*fs;
rx_fft = fft(rx_signals',Nfft)';
fft_mag = abs(rx_fft);
max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
[maxval,mindex] = max(fft_mag(:,max_search),[],2); % max in each row
carrier_phase = angle(rx_fft(max_search(mindex)'));
carrier_freq = fs/Nfft*max_search(mindex)';

carrier_freq = fc;

 % generate the time series and local oscillator
t = [0:1/fs:(rx_len-1)/fs];
lo = exp(1j*(2*pi*carrier_freq*t));
% downconvert
rx_baseband = rx_signals.*lo;
% slice out where data starts
rx_baseband = rx_baseband(1500:end);

% lowpass filtering both removes the 2fc term and anti-alias filters
% filtfilt used for 0 group delay filtering
rx_baseband = filtfilt(lpFilt,rx_baseband')';
expected_preamble1 = filtfilt(lpFilt,expected_preamble1')';
expected_preamble2 = filtfilt(lpFilt,expected_preamble2')';
expected_preamble_comb = filtfilt(lpFilt,expected_preamble_comb')';

% remove DC mean
rx_baseband = rx_baseband - mean(rx_baseband);
% rx_baseband = fftfilt(hpFilt,rx_baseband')';
% expected_preamble = filtfilt(hpFilt,expected_preamble);

sig_sec = rx_baseband;
Nfft = 2^nextpow2(length(sig_sec));
fft_sig = fft(sig_sec,Nfft);
[subcar_peak,subcar_peak_index] = max(fft_sig);
subcar_peak_f = subcar_peak_index/Nfft*fs;
subcar_peak_phase = angle(subcar_peak);
f = [-fs/2:fs/Nfft:fs/2-fs/Nfft];

t_window = 0.1;
window = chebwin(t_window*fs);
%Nfft = length(window);

% figure(1);
% hold on;
% %[pxx,f] = pwelch(sig_sec,window,[],Nfft,fs,'power');
% plot(f,20*log10(abs(fftshift(fft_sig))));
% grid on;
% grid minor;

% disp("Total Average Power (dB):");
% disp(num2str(10*log10(tot_pow)));

% figure(2);
% clf;
% hold on;
% plot(real(sig_sec(1:24000)));
% plot(imag(sig_sec(1:24000)));

ylim([-0.05 0.05]);

%% CORRELATE, ESTIMATE, DECODE, and COMPUTE BER %%

% N_data_bits = 0;
% 
% % resample for even integer fm0_samp value
% % with fs/fb divisble by 2 this code is not run
% if mod(fs/fb,2) ~= 0
%     p = (ceil(fs/fb)+1)*fb;
%     q = fs;
%     
%     rx_baseband = resample(rx_baseband,p,q);
%     
%     fs = p;
% end

fm0_samp = fs/fb;

% DO NOT CHANGE %
% data_len = round(N_data_bits/fb*fs);                % known data length in samples
preamble_len1 = round(length(expected_preamble1));    % known preamble length in samples
preamble_len2 = round(length(expected_preamble2));
packet_len = preamble_len1;                 % known packet length in samples

% number of packets to decode
% usually this needs to be set slightly lower than the total that was
% actually transmitted, I haven't investigated this further and am not sure
% why this is needed
N_packets1 = 4;
N_packets2 = 5;
N_comb_packets = 5;

% % estimates after channel projection for individual elements
% rx_rep_estimates = zeros(Nel,N_packets*(packet_len-preamble_len1));

channel_estimates = zeros(Nel,N_packets1+N_packets2+N_comb_packets);

A_hat = zeros(Nel,N_packets1+N_packets2+N_comb_packets);
ang_hat = zeros(Nel,N_packets1+N_packets2+N_comb_packets);

preamble_starts = zeros(Nel,1);

decode_preamble1 = expected_preamble1-mean(expected_preamble1);
decode_preamble2 = expected_preamble2-mean(expected_preamble2);
decode_preamble_comb = expected_preamble_comb-mean(expected_preamble_comb);

% figure;
% CORRELATION AND DECODING %
for pnum=1:N_packets1+N_packets2+N_comb_packets
    for el=1:Nel
        begdex = (pnum-1)*packet_len+1;
        endex = pnum*packet_len;
        end_preamble_dex = (pnum-1)*packet_len+1+preamble_len1;
    
%         beg_data_dex = (pnum-1)*data_len+1;
%         end_data_dex = pnum*data_len;
    
%         beg_bit_dex = (pnum-1)*N_data_bits+1;
%         end_bit_dex = pnum*N_data_bits;
        
        % perform cross correlation for packet start
        if pnum <= N_packets1
            decode_preamble = decode_preamble1;
        elseif pnum <= N_packets1+N_packets2
            decode_preamble = decode_preamble2;
        else
            decode_preamble = decode_preamble_comb;
        end
        
        % tx norm is length of preamble for binary keying (-1,+1)
        tx_norm = sum(abs(decode_preamble).^2);

        [rcorr,rlags] = xcorr(real(rx_baseband(el,begdex:end_preamble_dex))',decode_preamble');
        [icorr,ilags] = xcorr(imag(rx_baseband(el,begdex:end_preamble_dex))',decode_preamble');
        % removes tails of correlation 
        corr_tot = rcorr(end_preamble_dex-begdex+1:end)+1j*icorr(end_preamble_dex-begdex+1:end);
        abs_corr = abs(corr_tot);
        % find maximum correlation and begin decoding from there
        [preamble_max,preamble_starts(el)] = max(abs_corr); 

%         clf;
%         
%         plot(abs_corr/30);
%         hold on;
%         plot(real(rx_baseband(el,begdex:endex)));
%         plot([zeros(1,preamble_starts(el)) decode_preamble/10]);
        %xlim([0 1000]);

        % slice out packet found from correlation
        packet = rx_baseband(el,begdex+preamble_starts(el)-1:endex+preamble_starts(el)-1);
        % remove the mean
        packet = packet - mean(packet(1:preamble_len1));
        % slice out preamble
        packet_preamble = packet(1:preamble_len1);
        % slice out data
        packet_data = packet(preamble_len1+1:end);
        
        % compute the channel estimate 
        channel_estimates(el,pnum) = sum(packet_preamble.*conj(decode_preamble))/tx_norm;
        A_hat(el,pnum) = 2*abs(channel_estimates(el,pnum));
        ang_hat(el,pnum) = angle(channel_estimates(el,pnum));
        
        % extract estimate of data only
%         rx_rep_estimates(el,beg_data_dex:end_data_dex) = packet_data.*conj(channel_estimates(el,pnum))./abs(channel_estimates(el,pnum)).^2;
%         packet_preamble_proj = packet_preamble.*conj(channel_estimates(el,pnum));
%         preamble_estimate = packet_preamble_proj./abs(channel_estimates(el,pnum)).^2;
%         summed_preamble_estimate(pnum,:) = summed_preamble_estimate(pnum,:) + packet_preamble_proj;
%         % SNR calculation (currently incorrect)
%         noise_est = decode_preamble-preamble_estimate;
%         rx_snr_packet(el,pnum) = mean(abs(preamble_estimate).^2)/mean(abs(noise_est).^2);

%         plot(decode_preamble);
%         hold on;
%         plot(real(preamble_estimate));
        
%         % decode packet
%         if (mod(fm0_samp,2) == 0)
%             % should always enter this branch, camera_decode written by
%             % saad and waleed
%             bits = camera_decode(rx_rep_estimates(el,beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%         elseif mod(floor(fm0_samp),2)==0
%             bits = fm0_decode_new_R12_dec_even(rx_rep_estimates(el,beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%         else
%             bits = fm0_decode_new_R12(rx_rep_estimates(el,beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%         end
%         
%         % place bits into decoded data matrix
%         rx_decoded_data(el,beg_bit_dex:end_bit_dex) = bits;
    end
    
%     % sum array elements 
%     rx_rep_est_summed(beg_data_dex:end_data_dex) = sum(rx_rep_estimates(:,beg_data_dex:end_data_dex).*abs(channel_estimates(:,pnum)).^2,1)./sum(abs(channel_estimates(:,pnum)).^2,1);
%     summed_preamble_estimate(pnum,:) = summed_preamble_estimate(pnum,:) / sum(abs(channel_estimates(:,pnum)).^2,1);
%     summed_noise_est = decode_preamble-summed_preamble_estimate(pnum,:);
%     rx_snr_packet_summed(pnum) = mean(abs(summed_preamble_estimate(pnum,:)).^2)/mean(abs(summed_noise_est).^2);
% 
%     % decode summed packet
%     if (mod(fm0_samp,2) == 0)
%         bits = camera_decode(rx_rep_est_summed(beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%     elseif mod(floor(fm0_samp),2)==0
%         bits = fm0_decode_new_R12_dec_even(rx_rep_est_summed(beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%     else
%         bits = fm0_decode_new_R12(rx_rep_est_summed(beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%     end
%     
%     % place summed bits into decoded data matrix
%     rx_decoded_summed_data(beg_bit_dex:end_bit_dex) = bits;
end

mag_percent_error = abs(A_hat - Au*Ad) / (Au*Ad) * 100;
ang_percent_error = abs(ang_hat - angle(exp(1j*2*pi*fc*(tau_u1+tau_d1)))) / angle(exp(1j*2*pi*fc*(tau_u1+tau_d1))) * 100;

channel_mags_theta = [channel_mags_theta mean(A_hat(N_packets1+N_packets2+1:end))];

%% EXPORT DATA %%
% data_ch1 = data(t,preamble1,fb,n_data_reps1);
% data_ch2 = data(t-length(expected_preamble1)*n_data_reps1/fs-1e-3,preamble2,fb,n_data_reps2);
% 
% data_comb = data(t-(length(expected_preamble1)*n_data_reps1+length(expected_preamble2)*n_data_reps2)/fs-2e-3,preamble2,fb,n_data_reps2);
% 
% data_tot = (data_ch1+data_comb+1)/2.5+1j*(data_ch2+data_comb+1)/2.5;
% 
% figure;
% hold on;
% plot(real(data_tot));
% %plot(imag(data_tot));
% 
% write_complex_binary(data_tot,"../../tx_outputs/array_channel_estimating_data.dat");

theta = theta + 1;
channel_sim;

function out = data(t,code,fb,nreps)
    fm0_code = generate_fm0_sig2(code,2);
    code_len = length(fm0_code);
    out = zeros(size(t),'like',t);

    for n=1:nreps*code_len
        out(logical((t < n*1/(2*fb)).*(t >= (n-1)*1/(2*fb)))) = fm0_code(mod(n-1,code_len)+1);
    end
end
