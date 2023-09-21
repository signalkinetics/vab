imp1 = readmatrix("../../impedance/PAB010B_RX_IND+_1k-60k_801PTS_RIVER_ROTATOR_3MD_010B_008A.CSV");
imp1 = imp1(1:801,:);

imp2 = readmatrix("../../impedance/PAB008A_RX_IND+_1k-60k_801PTS_RIVER_ROTATOR_3MD_004A_008A.CSV");
imp2 = imp2(1:801,:);

fmin = 15e3;
fmax = 20e3;

f = linspace(fmin,fmax,1000);
w = 2*pi*f;

ls = 1.87e-3;

rlc1 = rlc_modeler(imp1,[fmin fmax],[],ls,1);
rlc2 = rlc_modeler(imp2,[fmin fmax],[],ls,1);

r1 = rlc1(1);
l1 = rlc1(2)*1e-3;
c1 = rlc1(3)*1e-9;
c01 = rlc1(4)*1e-9;
ls1 = ls;

r2 = rlc2(1);
l2 = rlc2(2)*1e-3;
c2 = rlc2(3)*1e-9;
c02 = rlc2(4)*1e-9;
ls2 = ls;

A1 = 1+1j*w*c01.*(1./(1j*w*c1)+1j*w*l1);
B1 = 1./(1j*w*c1)+1j*w*l1+1j*w*ls1+(1j*w*c01).*((1./(1j*w*c1)+1j*w*l1).*(1j*w*ls1));
C1 = 1j*w*c01;
D1 = 1+1j*w*ls1.*(1j*w*c01);

A2 = 1+1j*w*ls2.*(1j*w*c02);
B2 = 1./(1j*w*c2)+1j*w*l2+1j*w*ls2+(1j*w*c02).*((1./(1j*w*c2)+1j*w*l2).*(1j*w*ls2));
C2 = 1j*w*c02;
D2 = 1+1j*w*c02.*(1./(1j*w*c2)+1j*w*l2);

T1 = zeros(2,2,length(w));
T2 = zeros(2,2,length(w));

for n=1:length(w)
    T1(:,:,n) = [A1(n) B1(n);C1(n) D1(n)];
    T2(:,:,n) = [A2(n) B2(n);C2(n) D2(n)];
end

sij = get_sparams(T1,T2,[r1 r2]);

s11mag = squeeze(20*log10(abs(sij(1,1,:))));
s21mag = squeeze(20*log10(abs(sij(2,1,:))));

figure(3);
subplot(2,1,1);
hold on;
plot(f/1e3,s11mag);
title('S11 dB20');
xlabel("Freq (kHz)");
grid on;

subplot(2,1,2);
hold on;
plot(f/1e3,s21mag);
title('S21 dB20');
xlabel("Freq (kHz)");
grid on;
