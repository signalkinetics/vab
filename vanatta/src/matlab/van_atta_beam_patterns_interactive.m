global params phimin phimax ls;
%close all;
params.d = 7e-2;
params.c = 1500;
params.f = 18.5e3;
params.N = 4;
params.A = 1;
params.r = 1;
params.phi = 0;
params.do_direction_val = 0;

% imp1 = readmatrix("../../impedance/PAB008A_RX_IND+_1k-60k_801PTS_RIVER_ROTATOR_3MD_004A_008A.CSV");
% imp1 = imp1(1:801,:);
% 
% imp2 = readmatrix("../../impedance/PAB010B_RX_IND+_1k-60k_801PTS_RIVER_ROTATOR_3MD_010B_008A.CSV");
% imp2 = imp2(1:801,:);

dmin = 1e-2;
dmax = 30e-2;

phimin = -pi/2;
phimax = pi/2;

fmin = 10e3;
fmax = 30e3;
ls = 1.87e-3;
params.rlc1 = rlc_modeler(imp1,[fmin fmax],[],ls,1);
params.rlc2 = rlc_modeler(imp2,[fmin fmax],[],ls,1);

atot_db = generate_pattern(params);

global Ntheta theta min_r_db pol pol_steer;

min_r_db = 20*log10(params.A)-30;
max_r_db = 20*log10(params.A)+3+20*log10(params.N);
Ntheta = 1000;
theta = linspace(phimin,phimax,Ntheta); 

f = figure;
ax = axes('Parent',f,'position',[0.13 0.39 0.77 0.54]);
%subplot(1,3,3);
pol = polarplot(theta,atot_db);
hold on; 
pol_steer = polarplot(params.phi*ones(1,10),linspace(min_r_db,max_r_db,10),'--');
rlim([min_r_db max_r_db]);

global dslider_label phislider_label fslider_label;

%sliders and labels
dslider = uicontrol('Parent',f,'Style','slider','Position',[100,800,230,23],...
              'value',params.d, 'min',dmin, 'max',dmax,'Tag','d');
dslider_label = uicontrol('Parent',f,'Style','text','Position',[100,770,230,23],...
                'String',strcat('Element separation=',strcat(num2str(params.d/1e-2),'cm')));

phislider = uicontrol('Parent',f,'Style','slider','Position',[100,740,230,23],...
              'value',params.phi, 'min',phimin, 'max',phimax,'Tag','phi');
phislider_label = uicontrol('Parent',f,'Style','text','Position',[100,710,230,23],...
                'String',strcat('Angle=',strcat(num2str(params.phi*180/pi),'deg')));

fslider = uicontrol('Parent',f,'Style','slider','Position',[100,680,230,23],...
              'value',params.f, 'min',fmin, 'max',fmax,'Tag','f');
fslider_label = uicontrol('Parent',f,'Style','text','Position',[100,650,230,23],...
                'String',strcat('Frequency=',strcat(num2str(params.f/1e3),'kHz')));

Ndrop = uicontrol(f,'Style','popupmenu','Tag','N');
Ndrop.Position = [100,620,230,23];
Ndrop.String = {'2','4','6','8','10','12','14','16','18','20'};

Ndrop_label = uicontrol('Parent',f,'Style','text','Position',[100,590,230,23],...
                'String','N elements');

dslider.Callback = @update_plot; 
phislider.Callback = @update_plot; 
fslider.Callback = @update_plot; 
Ndrop.Callback = @update_plot;


% R1 = 400;
% L1 = 100e-6;
% C1 = 30e-9;
% 
% R2 = 65;
% L2 = 150e-6;
% C2 = 25e-9;
% 
% f = 20e3;
% s = 1j*2*pi*f;
% 
% y = 1/(1/(s*C1)+1/(s*C2)+s*(L1+L2));
% Y = [y -y; -y y];
% 
% F = [1/(2*sqrt(R1)) 0; 0 1/(2*sqrt(R2))];
% ZR = [R1 0; 0 R2];
% 
% i = eye(2);
% S = F*(i-ZR*Y)*(i+ZR*Y)^(-1)*F^(-1);
% 
% gamma1 = S(2,2);
% 
% Zr = R2;
% Zl = R1+1/(s*C1)+1/(s*C2)+s*(L1+L2);
% 
% gamma2 = (Zl-conj(Zr))/(Zl+Zr);

function atot = generate_pattern(params)
    global Ntheta theta min_r_db phimin phimax ls;

    at = zeros(params.N,Ntheta);    % transmitted wave column vector
    ar = zeros(params.N,Ntheta);    % reflected wave column vector
    an = zeros(params.N,Ntheta);    % total wave at observed point theta (far field)

    w = 2*pi*params.f;

    r1 = params.rlc1(1);
    l1 = params.rlc1(2)*1e-3;
    c1 = params.rlc1(3)*1e-9;
    c01 = params.rlc1(4)*1e-9;
    ls1 = ls;

    r2 = params.rlc2(1);
    l2 = params.rlc2(2)*1e-3;
    c2 = params.rlc2(3)*1e-9;
    c02 = params.rlc2(4)*1e-9;
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

    %sij = get_sparams(T1,T2,T1,T2,[r1 r2 r1 r2]);
    % MANUALLY ENTERED S PARAMS FROM ADS FOR 002A<->004D inner and
    % 001A<->004B outer POS STATE
%     sij = [ 0.798+1j*0.164 0 0 -0.12+1j*0.567;
%             0 0.748+1j*0.163 -0.155+1j*0.624 0;
%             0 -0.155+1j*0.624 0.738+1j*0.205 0;
%             -0.12+1j*0.567 0 0 0.796+1j*0.175];
    
    % 001A<->004B outer NEG STATE
%     sij = [ 0.798+1j*0.164 0 0 0.12-1j*0.567;
%             0 0.748+1j*0.163 0.155-1j*0.624 0;
%             0 0.155-1j*0.624 0.738+1j*0.205 0;
%             0.12-1j*0.567 0 0 0.796+1j*0.175];
    
     sij = [ 0 0 0 1;
            0 0 1 0;
            0 1 0 0;
            1 0 0 0];
%     sij = [0 1;
%            1 0];

    lambda = params.c / params.f;
    k = 2*pi/lambda;

    phi = params.phi;

    if params.do_direction_val == 1
         phi = theta;
    end
    
    for n=1:params.N
        rn = exp(1j*k*(params.r-(n-1)*params.d*sin(theta)));
        at(n,:) = params.A*(sij(n,params.N+1-n)*exp(-1j*k*sin(phi)*(params.N-n)*params.d)).*rn;
        ar(n,:) = params.A*(sij(n,n)*exp(-1j*k*sin(phi)*(n-1)*params.d)).*rn;
    end
    an = ar+at;
    atot = abs(sum(an,1));
    %atot_norm = atot/max(atot);
    atot = 20*log10(atot);

    phase_deg = angle(sum(an,1))*180/pi;

    direction_val = atot(floor((params.phi-phimin)/(phimax-phimin)*(Ntheta-1)+1));
    phase_direction_val = phase_deg(floor((params.phi-phimin)/(phimax-phimin)*(Ntheta-1)+1));

    disp(['Value in direction: ' num2str(direction_val) 'dB']);
    disp(['Phase in direction: ' num2str(phase_direction_val) ' deg']);

    atot(atot < min_r_db) = min_r_db;
end

function a = update_plot(es,ed)
    global dslider_label phislider_label fslider_label pol pol_steer params;

    switch ed.Source.Tag
        case 'd'
            params.d = es.Value;
            dslider_label.String = strcat('Element separation=',strcat(num2str(es.Value/1e-2),'cm'));
        case 'phi'
            params.phi = es.Value;
            phislider_label.String = strcat('Angle=',strcat(num2str(es.Value*180/pi),'deg'));
        case 'f'
            params.f = es.Value;
            fslider_label.String = strcat('Frequency=',strcat(num2str(es.Value/1e3),'kHz'));
        case 'N'
            params.N = 2*es.Value;
        otherwise
    end

    set(pol,'YData',generate_pattern(params));
    set(pol_steer,'XData',params.phi*ones(1,10));
end


