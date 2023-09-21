fs = 2e5;
fc=20;

amp = 0.3;

t_length = 1;
t_pulse = 10e-3;
t = [0:1/fs:t_length-1/fs];

y = amp*(t < t_pulse);

write_complex_binary(y, strrep("../../tx_outputs/pulse_?sec.dat","?",num2str(t_length)));
