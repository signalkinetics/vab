fs = 2e5;

root = '../../rx_outputs/WHOI Van Atta 2 Microbenchmarks 12-01-2022/';
filename = 'fixed_array_chest_006A_006C_txfmr_nicktb_siggen_18,5kfc_0,0deg_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_3m_depth_005B_purui_tx_60Vrms_1,9m_1m_hphydro_diff_0.dat';
filepath = strcat(root,filename);

yr = read_complex_binary(filepath);        
sig = yr(24:end);
sig = real(sig)-imag(sig);

N = length(sig);
t_tot = N/fs;

fc = 18.5e3;

ph_err = [];

t = zeros(1,N);
vco = zeros(1,N);
ph = zeros(1,N);
ph_est = zeros(1,N);
lp = zeros(1,N);
y = zeros(1,N);
integ = zeros(1,N);

Bn = 0.01*fs;
damp = 1;

k0 = 1;
kd = 1;
kp = 1/(kd*k0)*4*damp/(damp+1/(4*damp))*Bn/fs;
ki = 1/(kd*k0)*4/(damp+1/(4*damp))^2*(Bn/fs)^2;

integ_out = 0;

for n = 1:N-1
    t(n) = t_tot*n/N;
    % input signal
    y(n) = sig(n);%sin(2*pi*fc*t(n)+pi);

    % phase detect
    ph(n) = kd*y(n)*imag(vco(n));

    % loop filter
    integ_out = ki*ph(n)+integ_out;
    lp(n) = kp*ph(n) + integ_out;

    % vco
    ph_est(n+1) = ph_est(n) + k0*lp(n);
    vco(n+1) = exp(-1j*(2*pi*fc*t_tot*(n+1)/N+ph_est(n)));
end


plot(t(1:end-1),y(1:end-1));
hold on;
plot(t(1:end-1),real(vco(1:end-1)));