fs = 2e5;
fc=18.5e3;

amp = 230e-3;

t_length = 35;
t = [0:1/fs:t_length-1/fs];

y = amp*sin(2*pi*fc*t);

write_complex_binary(y, strrep("../../tx_outputs/tone_18,5kfc_460mVpp_?sec.dat","?",num2str(t_length)));
