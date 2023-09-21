folder = '~/Documents/sk/oceans/vanatta/rx_outputs/River PAB Phase Tests 06-08-2022/';
file = 'rx_phase_pab_007B_006A_007A_004A_*deg_2m_depth_3,5m_u2b';
root = strcat(folder,file);

node_list = ["007B","006A","007A","004A"];

fs = 2e5;
Nsamps = 10*fs;

bp = designfilt('bandpassfir', ...       % Response type
       'StopbandFrequency1',17e3, ...    % Frequency constraints
       'PassbandFrequency1',18e3, ...
       'PassbandFrequency2',19e3, ...
       'StopbandFrequency2',20e3, ...
       'DesignMethod','kaiserwin', ...         % Design method
       'StopbandAttenuation1',40, ...         % Design method options
       'StopbandAttenuation2',40, ...
       'PassbandRipple',1, ...
       'SampleRate',fs);               % Sample rate

d = 0.062;
lambda = 1500/18.5e3;

t_window = 10;
window = chebwin(t_window*fs);
Nfft = length(window);
i = 1;

for ang=-45:45:45
    exp_phase_diff = [0];
    act_phase_diff = [0];

    filename = strrep(root,"*",num2str(ang));

    data = read_complex_binary(strcat(filename,"_0",".dat"));
    node1 = real(data(1:Nsamps));
    node1 = filter(bp,node1);
    prev_meas_ang = 0;

    for n=2:length(node_list)
        read_filename = strcat(filename,"_",num2str(floor((n-1)/2)),".dat");
        data = read_complex_binary(read_filename);
        if mod(n,2) == 0
            node2 = imag(data(1:Nsamps));
        else
            node2 = real(data(1:Nsamps));
        end
   
        node2 = filter(bp,node2);
        
        fft1 = fft(node1.*window,Nfft);
        fft2 = fft(node2.*window,Nfft);
        
        fft1 = fft1(1:Nfft/2);
        fft2 = fft2(1:Nfft/2);
        
        max1 = max(fft1);
        max2 = max(fft2);
        
        meas_ang = angle(max2/max1);

        if ang < 0
            while meas_ang > prev_meas_ang
                meas_ang = meas_ang - 2*pi;
            end
        end
        if ang > 0
            while meas_ang < prev_meas_ang
                meas_ang = meas_ang + 2*pi;
            end
        end

%         if n >=2 && ang ~= 0
%             meas_ang = meas_ang + 2*pi*(n-1)*sign(ang);
%         end

        exp_phase_diff = [exp_phase_diff 2*pi*(n-1)*d*sin((ang-4)/180*pi)/lambda];
        act_phase_diff = [act_phase_diff meas_ang];

        prev_meas_ang = meas_ang;
    end

    subplot(3,1,i);
    plot(exp_phase_diff/pi*180);
    hold on;
    plot(act_phase_diff/pi*180);
    xlabel("n");
    ylabel("Phase (deg)");
    title(strcat("Phase vs. element (no offset) at ",num2str(ang)," deg"));
    legend("Expected","Actual");
    %ylim([-180 180]);
    xlim([1 length(node_list)]);
    grid on;

    i = i+1;
end

disp("Expected phase difference (deg): ");
disp(num2str(exp_phase_diff/pi*180));

disp("Actual phase difference (deg): ");
disp(num2str(act_phase_diff/pi*180));