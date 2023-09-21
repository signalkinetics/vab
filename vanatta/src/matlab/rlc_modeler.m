function rlc_c0 = rlc_modeler(measured_data,fminmax,init_rlc_c0,ls,display)

    fmin = fminmax(1);
    fmax = fminmax(2);
    
    f = measured_data(:,1);
    
    ifmin = find(f >= fmin,1);
    ifmax = find(f <= fmax,1,'last');
    
    f = f(ifmin:ifmax);
    r = measured_data(ifmin:ifmax,2);
    x = measured_data(ifmin:ifmax,3);
    ydata = [r';x'];
    zdata = r+1j*x;
    % % % 
%     n = 2;
    if isempty(init_rlc_c0)
        rlc0(1) = 300; % ohm
        rlc0(2) = 6; % mH
        rlc0(3) = 15; % nF
        rlc0(4) = 40; % nF

        rlc0 = rlc0';
    else
        rlc0 = init_rlc_c0;
    end
    
    func_min = @(rlc) rlc_model(rlc,f)-ydata;
    
    options = optimoptions('lsqnonlin','FunctionTolerance',1e-9,'StepTolerance',1e-9,'MaxFunctionEvaluations',10e3,'Display','iter');
    [rlc,~,~,~,~] = lsqnonlin(func_min,rlc0,[0 0 0 0]',[],options);
    
    %rlc = rlc0;
    imp_opt = rlc_model(rlc,f);
    zopt = imp_opt(1,:)+1j*imp_opt(2,:);
    
    if display
        figure;
        subplot(1,2,1);
        plot(f/1e3,real(zdata));
        hold on;
        plot(f/1e3,real(zopt));
        title("Measured Real Part vs. LSQ RLC Real Part Fit");
        legend("Measured","LSQ Fit");
        xlabel("Frequency (kHz)");
        ylabel("Real Part (Ohm)");
        
        subplot(1,2,2);
        plot(f/1e3,imag(zdata));
        hold on;
        plot(f/1e3,imag(zopt));
        title("Measured Imaginary Part vs. LSQ RLC Imaginary Part Fit");
        legend("Measured","LSQ Fit");
        xlabel("Frequency (kHz)");
        ylabel("Imaginary Part (Ohm)");
    
        disp(strcat("R = ",num2str(rlc(1)), " ohm"));
        disp(strcat("L = ",num2str(rlc(2))," mH"));
        disp(strcat("C = ",num2str(rlc(3))," nF"));
        disp(strcat("C0 = ",num2str(rlc(4))," nF"));
        %disp(strcat("Ls = ",num2str(rlc(5))," mH"));
    end

    rlc_c0 = rlc;

    function yout = rlc_model(rlc,fdata)

        wdata = 2*pi*fdata;
        yout = zeros(2,length(fdata)); % allocate yout
    
        R1 = rlc(1);
        L1 = rlc(2)*1e-3;
        C1 = rlc(3)*1e-9;
        C0 = rlc(4)*1e-9;
        %L2 = rlc(5)*1e-3;
        
    %     yout(:,1) = real(R1+1j*wdata*L1+1./(1j*wdata*C1));
    %     yout(:,2) = imag(R1+1j*wdata*L1+1./(1j*wdata*C1));
        imp = 1j*wdata*ls+(R1+1j*wdata*L1+1./(1j*wdata*C1))./(1+1j*wdata*C0.*(R1+1j*wdata*L1+1./(1j*wdata*C1)));

        yout(1,:) = real(imp);
        yout(2,:) = imag(imp);
    end
end
