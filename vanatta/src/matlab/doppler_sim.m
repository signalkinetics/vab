fs = 2e5;
fb = 500;
fm0_samp = fs/fb;

fc = 18.5e3;
fc_h = 18.5e3;
fd = 0.05;
dmax = 0.01;
d = 10;

c = 1500;

Ab = 0.01;

lambda = c/fc;

preamble = [0 0 1 1 1 0 1 0];
%preamble = [0 0 1 1 1 0 1 0 0 1];
%preamble = [0 0 1 1 1 0 1 0 0 1 0 0 0 1 1 1];
expected_preamble = generate_fm0_sig2(preamble,fm0_samp);
    
N_preamble_bits = length(preamble);
N_data_bits = 16;

N_packets = 625;
packet_delay = 0; % delay in between each packet
N_tot_data_bits = N_data_bits*N_packets;
N_tot_bits = (N_data_bits+N_preamble_bits)*N_packets;

expected_data = real(read_complex_binary(strrep('../../../tx_outputs/data_prbs_order=15_len=?_packets=625.dat','?',num2str(N_data_bits))))';

% %%%%% RESHAPE EXPECTED DATA INTO FORMAT FOR DFE %%%%%%
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
t = [0:1/fs:(length(expected_data_signal)-1)/fs];

dd = dmax*sin(2*pi*fd*t);
dv = 2*pi*fd*dmax*cos(2*pi*fd*t);

re_xr = 1/2*(cos(2*pi*(fc_h-fc)*t)+Ab*expected_data_signal.*cos(2*pi*(fc_h-fc*(1-dv/c).^2).*t+4*pi*(d+dd)/lambda));

plot(t,re_xr);