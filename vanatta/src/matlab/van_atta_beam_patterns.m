
d = 3e-2;
c = 1500;
f = 18.5e3;
N = 2;
A = 1;
r = 1;
phi = pi/6;
do_direction_val = 1;

% imp1 = readmatrix("../../impedance/PAB008A_RX_IND+_1k-60k_801PTS_RIVER_ROTATOR_3MD_004A_008A.CSV");
% imp1 = imp1(1:801,:);
imp2 = 40-1j*10;


% imp2 = readmatrix("../../impedance/PAB010B_RX_IND+_1k-60k_801PTS_RIVER_ROTATOR_3MD_010B_008A.CSV");
% imp2 = imp2(1:801,:);
imp2 = 40+1j*10;

dmin = 1e-2;
dmax = 30e-2;

phimin = -pi/2;
phimax = pi/2;

fmin = 18.5e3;
fmax = 18.5e3;

imp1 = [fmin 40 10];
imp2 = [fmin 40 10];
% ls = 1.87e-3;
ls = 0;
rlc1 = rlc_modeler(imp1,[fmin fmax],[],ls,1);
rlc2 = rlc_modeler(imp2,[fmin fmax],[],ls,1);

% atot_db = generate_pattern(params);

min_r_db = 20*log10(A)-30;
max_r_db = 20*log10(A)+3+20*log10(N);
Ntheta = 1000;
theta = linspace(phimin,phimax,Ntheta); 


at = zeros(N,Ntheta);    % transmitted wave column vector
ar = zeros(N,Ntheta);    % reflected wave column vector
an = zeros(N,Ntheta);    % total wave at observed point theta (far field)

w = 2*pi*f;

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
B1 = 1./(1j*w*c1)+1j*w*l1+1j*w*ls1+1j*w*c01.*((1./(1j*w*c1)+1j*w*l1).*1j*w*ls1);
C1 = 1j*w*c01;
D1 = 1+1j*w*ls1.*(1j*w*c01);

A2 = 1+1j*w*ls2.*(1j*w*c02);
B2 = 1./(1j*w*c2)+1j*w*l2+1j*w*ls2+1j*w*c02.*((1./(1j*w*c2)+1j*w*l2).*1j*w*ls2);
C2 = 1j*w*c02;
D2 = 1+1j*w*c02.*(1./(1j*w*c2)+1j*w*l2);

T1 = [A1 B1;C1 D1];
T2 = [A2 B2;C2 D2];

sij = get_sparams(T1,T2,[r1 r2]);
% MANUALLY ENTERED S PARAMS FROM ADS FOR 002A<->004D inner and
% 001A<->004B outer POS STATE
%     sij = [ 0.798+1j*0.164 0 0 -0.12+1j*0.567;
%             0 0.748+1j*0.163 -0.155+1j*0.624 0;
%             0 -0.155+1j*0.624 0.738+1j*0.205 0;
%             -0.12+1j*0.567 0 0 0.796+1j*0.175];

%     sij = [ 0.798+1j*0.164 0 0 0.12-1j*0.567;
%             0 0.748+1j*0.163 0.155-1j*0.624 0;
%             0 0.155-1j*0.624 0.738+1j*0.205 0;
%             0.12-1j*0.567 0 0 0.796+1j*0.175];
%sij = [sqrt(2)/2 sqrt(2)/2; sqrt(2)/2 sqrt(2)/2];

lambda = c / f;
k = 2*pi/lambda;

if do_direction_val == 1
     phi = theta;
end

for n=1:N
    rn = exp(1j*k*(r-(n-1)*d*sin(theta)));
    at(n,:) = A*(sij(n,N+1-n)*exp(-1j*k*sin(phi)*(N-n)*d)).*rn;
    ar(n,:) = A*(sij(n,n)*exp(-1j*k*sin(phi)*(n-1)*d)).*rn;
end
an = ar+at;
atot = abs(sum(an,1));
%atot_norm = atot/max(atot);
atot_db = 20*log10(atot);

at_db = 20*log10(abs(sum(at,1)));
ar_db = 20*log10(abs(sum(ar,1)));

phase_deg = angle(sum(an,1))*180/pi;

direction_val = atot_db(floor((phi-phimin)/(phimax-phimin)*(Ntheta-1)+1));
phase_direction_val = phase_deg(floor((phi-phimin)/(phimax-phimin)*(Ntheta-1)+1));

disp(['Value in direction: ' num2str(direction_val) ' dB']);
disp(['Phase in direction: ' num2str(phase_direction_val) ' deg']);

atot_db(atot_db < min_r_db) = min_r_db;
at_db(at_db < min_r_db) = min_r_db;
ar_db(ar_db < min_r_db) = min_r_db;

f = figure;
%pol = polarplot(theta,atot_db);
% hold on; 

polarplot(theta,at_db, 'LineWidth',6);
hold on;
polarplot(theta,ar_db,'--','LineWidth',6);
% pol_steer = polarplot(phi*ones(1,10),linspace(min_r_db,max_r_db,10),'--','Color','blue');
% pol_steer = polarplot(-phi*ones(1,10),linspace(min_r_db,max_r_db,10),'--','Color','red');
rlim([min_r_db max_r_db]);
ax = gca;
% ax.FontSize = 36;
% rule = ax.RAxis;
% rlabel = rule.Label;
% rlabel.FontSize = 36;
% rlabel.Rotation = 90;
rlabel.String = "Magnitude (dB) \rightarrow";
legend("Transmission Response","Reflection Response");



