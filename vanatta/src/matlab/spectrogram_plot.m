folder = '~/Documents/sk/oceans/vanatta/rx_outputs/River PAB Van Atta 06-01-2022/';
file = 'rx_river_backscatter_pab_007A_007B_ind_vanatta_-45deg_tmux_18,5kfc_1kHz_square_1m_depth_5,8m_dis_4,8m_hphydro_120sec_0.dat';

root = strcat(folder,file);

fs = 2e5;

sig = real(read_complex_binary(root));
%sig = sig(1:fs*120);
%sig = cos(2*pi*20e3*linspace(0,120,fs*120));

t_window = 1;
window = chebwin(t_window*fs);
Nfft = length(window);

[pxx,f] = pwelch(sig,window,[],Nfft,fs,'power');
figure;
plot(f(17e3*Nfft/fs:20e3*Nfft/fs),10*log10(pxx(17e3*Nfft/fs:20e3*Nfft/fs)));

%peak_bin_pow = max(pxx)*fs/Nfft_welch;

% 
% Nfft = length(sig(1:fs*10));
% ideal_fft = fft(sig(1:fs*10),Nfft)/Nfft;
% f_fft = 0:fs/Nfft:fs/2-fs/Nfft;
% 
% hold on;
% plot(f_fft,20*log10(abs(ideal_fft(1:Nfft/2))));



