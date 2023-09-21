Npts = 801;

node1 = readmatrix("../../impedance/PAB011A_RX_1k-60k_801PTS_2,5MD_RIVER_ROTATOR_008A_011B_011A_010B.CSV");
node1 = node1(1:Npts,:);

node2 = readmatrix("../../impedance/PAB010B_RX_1k-60k_801PTS_2,5MD_RIVER_ROTATOR_008A_011B_011A_010B.CSV");
node2 = node2(1:Npts,:);

% node3 = readmatrix("../../impedance/PAB4_RX_1k-60k_IMP_RIVER_ROTATOR_8020.CSV");
% node3 = node3(1:201,:);
% 
% node4 = readmatrix("../../impedance/PAB4G_IND_RX_IMPEDANCE_SG.CSV");
% node4 = node4(1:201,:);

N=2;

% find zero crossings
zerox = zeros(N,Npts);
real_zerox = zeros(N,Npts);

figure(1);
lgd1 = legend;

figure(2);
lgd2 = legend;

for n=1:N
    switch n
        case 1
            node = node1;
        case 2
            node = node2;
        case 3
            node = node3;
        case 4
            node = node4;
        otherwise
            node = node1;
    end

    for i=2:length(node1)
        if node(i,3)*node(i-1,3) < 0
            x1 = node(i-1,1);
            y1 = node(i-1,3);
            x2 = node(i,1);
            y2 = node(i,3);

            y1_r = node(i-1,2);
            y2_r = node(i,2);

            zerox(n,i) = x1-y1*(x2-x1)/(y2-y1); % calculates proper x value (frequency) at which linearly interpolated ZC occurs
            real_zerox(n,i) = y1_r+(y2_r-y1_r)/(x2-x1)*(zerox(n,i)-x1);
        end
    end

    figure(1);
    hold on;
    plot(node(:,1)/1e3,node(:,2));
    lgd1.String{n} = strcat('PAB ',num2str(2*n),' R');
    xlabel("Frequency (kHz)");
    ylabel("Real Part (Ohm)");
    
    zerox_plot = nonzeros(zerox(n,:))/1e3;

    figure(2);
    hold on;
    p = plot(node(:,1)/1e3,node(:,3));
    text(zerox_plot,zeros(1,length(zerox_plot)),'x','Color',p.Color);

    lgd2.String{n} = strcat('PAB ',num2str(2*n),' X');
    xlabel("Frequency (kHz)");
    ylabel("Imaginary Part (Ohm)");
    
end

figure(1);
grid on;

figure(2);
grid on;

% legend('PAB2 R','PAB2 X');
% xlabel("Frequency (kHz)");
% ylabel("Value (Ohm)");