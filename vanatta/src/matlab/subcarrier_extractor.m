clear pxx;

% place question mark where degree # should be
folder = '~/Documents/sk/oceans/vanatta/rx_outputs/River PAB Van Atta 06-23-2022/';
file = "rx_backscatter_array_007B_003A_purui_match_0deg_tmux_21,2kfc_1kmod_2m_depth_2m_u2b_1m_midpower_hphydro_diff";
root = strcat(folder,file);

fmod = 1e3;
fs = 2e5;
Nsamps = 10*fs;

degree_list = [0];
Ndeg = length(degree_list);

measured_pattern = zeros(Ndeg,2);
measured_pattern(:,1) = degree_list * pi/180;

trial_length = fs*0.1;
Ntrials = 2*floor(Nsamps/trial_length);

avg_subcarrier_pow = zeros(Ntrials-1,Ndeg);
rx_signals = zeros(Nsamps,Ndeg);

for n=1:Ndeg
    filename = strrep(root,"?",num2str(degree_list(n)));
    
    sig = read_complex_binary(strcat(filename,'_0','.dat'));
    rx_signals(:,n) = real(sig(1:Nsamps));
end

for i=1:Ntrials-1
    segment = rx_signals((i-1)*trial_length/2+1:(i+1)*trial_length/2,:);
    Nfft = trial_length;
    [pxx,f] = periodogram(segment,chebwin(trial_length),Nfft,fs,'power');
    %pxx = pxx*fs/Nfft;
%     Nfft = trial_length;
%     seg_fft = fft(segment,Nfft);
%     pxx = abs(seg_fft(1:Nfft/2,:)).^2;
%     f = 0:fs/Nfft:fs/2-fs/Nfft;
    
    [carrier_pow,carrier_index] = max(pxx);
    subcarrier_left_index = carrier_index - fmod*Nfft/fs;
    subcarrier_right_index = carrier_index + fmod*Nfft/fs;
    
    subcarrier_left_window = linspace(subcarrier_left_index(1) - fmod/4*Nfft/fs,subcarrier_left_index(1) + fmod/4*Nfft/fs,fmod/2*Nfft/fs+1);
    subcarrier_right_window = linspace(subcarrier_right_index(1) - fmod/4*Nfft/fs,subcarrier_right_index(1) + fmod/4*Nfft/fs,fmod/2*Nfft/fs+1);
    
    [subcarrier_left_max,subcarrier_left_mindex] = max(pxx(subcarrier_left_window,:),[],1);
    [subcarrier_right_max,subcarrier_right_mindex] = max(pxx(subcarrier_right_window,:),[],1);
    
    avg_subcarrier_pow(i,:) = mean(cat(1,subcarrier_left_max,subcarrier_right_max));

    if 10*log10(avg_subcarrier_pow(i,1)) > -85
        test = 0;
    end

%     plot(f,pxx);
end

figure(1);
for n=1:Ndeg
    [cdf,x] = ecdf(10*log10(avg_subcarrier_pow(:,n)));
    hold on;
    plot(x,cdf);
    %title(strcat("ECDF of ",num2str(degree_list(n)),' deg'));
    grid on;

    disp(strcat("Avg Pow @ ",num2str(degree_list(n))," deg (dB): "));
    disp(num2str(mean(10*log10(avg_subcarrier_pow(:,n)))));
end

legend(num2str(degree_list'));

% xlim([-120 -80]);

%     if Ntrials > 1
%         for fnum=0:Ntrials-1
%             sig = read_complex_binary(strcat(filename,'_',int2str(fnum),'.dat'));
%             
%             rx_signals(:,fnum+1) = real(sig);
%         end
%     else
%         sig = read_complex_binary(strcat(filename,'_0.dat'));
%         rx_signals(:,1) = real(sig);
%     end
%     
   
%     avg_subcar_pow_trial(n,:) = mean(cat(1,subcarrier_left_max,subcarrier_right_max),1);
    
%     measured_pattern(n,2) = 10*log10(avg_subcarrier_pow);
%     disp("Average subcarrier power (dB): ");
%     disp(10*log10(avg_subcarrier_pow));
    



% figure(1);
% polarplot(measured_pattern(:,1),measured_pattern(:,2));
% hold on;
% rlim([-80 -50]);
% 
% figure(2);
% xlabel("Subcarrier Power (dB)");
% ylabel("F(x)");



% plot(f,10*log10(pxx(:,1)));

