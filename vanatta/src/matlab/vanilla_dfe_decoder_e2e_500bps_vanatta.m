clear all;
clc;
Fs = 2e5;
%path = '/Users/nazishnaeem/uhd/build/examples/';
path  = ['../../rx_outputs/River PAB Van Atta 4 09-16-2022/'];
%filename = 'rx_DFE_test_100m_500bps0.dat';
%filename = 'rx_h_3m_b_3m_t_500_6July_3_4_power0.dat';
filename = 'rx_single_chest_pab_012A_stag9cm_7cm_sp_2,9mtxfmr_+0deg_mosfet_18,5kfc_prbs_0,5kbps_usrp_2,5m_depth_010A_purui_new_tx_6m_5m_hphydro_400mVpp_0.dat';
% filenghahgme = 'test0.dat';
data2 = read_complex_binary([path filename]);


%%
data2 = [data2;zeros(300*Fs,1)];
Y = fftshift((fft(data2)));
F = -Fs/2:Fs/length(Y):Fs/2 - Fs/length(Y);
data_rate = 500; %10e3
fm0_samp = Fs/data_rate;
%fm0_samp = 195;
[val,args] = max(abs(Y));
freq = abs(F(args)); %Frequency at which the signal is centered
%freq = 20e3;
ph_c = angle(Y(args));

disp(['Carr Frequency: ' num2str(freq)]) 
data = data2(1:end-300*Fs);
Y2 = fftshift((fft(data)));
F2 = -Fs/2:Fs/length(Y2):Fs/2 - Fs/length(Y2);
[val2,args2] = max(abs(Y2));
freq2 = abs(F2(args2));
lpb = 4000; 
lsb = 5000;
b_cutoff = data_rate*5;

if b_cutoff > freq
    b_cutoff = 18000;
end

t = (0:1/Fs:length(data)/Fs - 1/Fs); % creating time vector
mult_s = ((data.').*exp(-1j*(2*pi*freq.*t-ph_c)));


FIRLPF = dsp.LowpassFilter('SampleRate',Fs,'PassbandFrequency',lpb, 'StopbandFrequency',lsb);
%FIRLPF = dsp.LowpassFilter('SampleRate',Fs,'PassbandFrequency',2300, 'StopbandFrequency',5000);

shifted_sig = FIRLPF([1 zeros(1,Fs)].');
[max_filt,del_indxlp] = max(shifted_sig);
disp(['Low Pass Filter Delay(Samples): ' num2str(del_indxlp)]);
release(FIRLPF)


final_lpf = [];
chunk = 1;
mult_s = mult_s(1:length(mult_s) - mod(length(mult_s),chunk));
final_lpf = [];

for i = 1:length(mult_s)/chunk:length(mult_s)

    pad = zeros(1*length(mult_s)/chunk,1);
    temp_lpf = FIRLPF([pad;mult_s(i:i+length(mult_s)/chunk-1).';pad]);
    temp2_lpf = temp_lpf(length(pad)+del_indxlp:end-length(pad)+del_indxlp-1);
    final_lpf = [final_lpf;temp2_lpf];
    %i
    
end

x = final_lpf;

xx = x;
Fpass = 150; %220;%1200;%80;
Fstop = 40; %20;%1000;%20;
filtertype = 'FIR';
% Astop = 80;
FIRHPF = dsp.HighpassFilter('SampleRate',Fs,...
                              'FilterType',filtertype,...
                              'PassbandFrequency',Fpass,...
                              'StopbandFrequency',Fstop);
                          
shifted_sig = FIRHPF([1 zeros(1,Fs)].');
[max_filt,del_indx] = max(shifted_sig);
disp(['Filter Delay(Samples): ' num2str(del_indx)]);
release(FIRHPF)
final_hpf = [];
chunk = 1;
xx = xx(1:length(xx) - mod(length(xx),chunk));
for i = 1:length(xx)/chunk:length(xx)
   
    pad = zeros(length(xx)/chunk,1);
    temp_hpf = FIRHPF([pad;xx(i:i+length(xx)/chunk-1);pad]);
    temp2_hpf = temp_hpf(length(pad)+del_indx:end-length(pad)+del_indx-1);
    final_hpf = [final_hpf;temp2_hpf];
    i
    
end
% 
final_hpf1 = my_lpf(final_hpf.', lpb,lsb,Fs);

final_hpf = final_hpf1;
final_hpf1 = my_hpf(final_hpf,Fpass,Fstop,Fs);
final_hpf = final_hpf1;

%%
%load bit_data
load jack_data_vanatta;
%bits = [complete_bits(1,:),complete_bits(2,:),complete_bits(3,:) complete_bits(4,:)];
pt = 20;
bits = [];
for q = 1:pt
  bits = [bits, complete_bits(q,:)];
end

    
pkt_size = 50;

fm0_samp = Fs/data_rate; %96.6; %96.6; %100.6; %190;
%load data_rate1000;
first_bits = bits(1:10);
first_pre = generate_fm0_sig2(first_bits,fm0_samp);
first_pre = generate_fm0_sig2(bits,fm0_samp);
mid_pre = first_pre;

final_hpf2 = final_hpf(1:(pt+20)*pkt_size*fm0_samp + 15e5);
%final_hpf2 = final_hpf(1:15e5);
%final_hpf2 = res(1:1e5);
% 
zz_r1 =  normxcorr2(first_pre.',real(final_hpf2));
zz_r22 = normxcorr2(first_pre.',imag(final_hpf2));
% 
zz_r1 = zz_r1(length(first_pre):end-length(first_pre)+1);
zz_r22 = zz_r22(length(first_pre):end-length(first_pre)+1);
% 
z_tot1 = abs(zz_r1+1i*zz_r22);  
[val_max,ind] = max(z_tot1);
%ind = inddd.Position(1);
plot(z_tot1)
y_rx = final_hpf(ind:end-0.5*Fs);

disp("Before DFE");
% 
% save('dfe_ofdm_sim_carr_0_045','y_rx')
% 
%%

%clc;
%close all;
%clear all;


%load dfe_ofdm_river_5kbps_3m
%load dfe_ofdm_river_5kbps_3m_real_tightfilt


%load dfe_ofdm_river_5kbps_10m_real_chunk1
%load rx_ofdm_chirp_5m_2phatay

clearvars -except y_rx weights
num_pkt_ind =100;
num_testing_pkts = 350;
Fs = 2e5;
W = 0;
%W = weights;
load jack_data_vanatta;
ds_ind = 1;%1;%5;
y_rx = downsample(y_rx,ds_ind);
y_rx = y_rx(1:end-0.1*Fs);
f_factor = 0.999;%[0.98 0.983 0.989 0.999 1];
%num_pkt_ind = 1;%[90:1:105];
ref_tap = 1;ff_tap = 203;%53;%23; 
fb_tap =100;%200;%100;
const_lvl =1.126e-4;%max(real(y_rx(3.84e6:3.89e6)));%(1.0350e-04); %(3.126e-4);%3.126e-4; 
%num_pkt_ind = 28*ones(1,20);
Fs = 2e5/ds_ind;
data_rate = 500;%5000;%2e3;
fm0samp = Fs/data_rate;
pkt_bits = 24; %50 for waleed 24 for jack
pre_bits = 8; % 10 for waleed 8 for jack
pkt_size = pkt_bits*fm0samp;



row = 1;
for ff = f_factor
    %clearvars -except ff num_pkt row num_pkt_ind f_factor y_rx complete_bits ber_vec_a ber_vec_b cnttt row 
    
    chunk = 1;
    %num_pkt = 5;
    off = 1;
    temp = chunk;
    training_pkt = [];
    training_bits = [];
    
    cnttt = 1;
    
    for num_pkt = num_pkt_ind
        training_pkt = [];
        training_bits = [];
        temp = chunk;
        for ss = 1:num_pkt
%             if(num_testing_pkts == 26)
%                 disp('ola')
%             end
            training_bits = [training_bits complete_bits(temp,:)];
            training_pkt = [training_pkt generate_fm0_sig2([complete_bits(temp,:)],fm0samp).*const_lvl]; 
            temp = temp + 1;
        end
        
        train_symb = downsample(training_pkt,fm0samp/2);

        rx_data = [y_rx((chunk-off)*pkt_size+1:(chunk-off)*pkt_size+num_pkt*pkt_size+num_testing_pkts*pkt_size)];

         acorr = 0.1;
         eqdfe_lms = comm.DecisionFeedbackEqualizer('Algorithm','RLS','InputSamplesPerSymbol',fm0samp/2, 'ReferenceTap',ref_tap, ...
             'ForgettingFactor',ff,'NumForwardTaps',ff_tap,'AdaptAfterTraining',1,'NumFeedbackTaps',fb_tap,...
             'Constellation',[min(training_pkt) max(training_pkt)],'InitialInverseCorrelationMatrix',acorr,'InitialWeights',(W),'InitialWeightsSource','Property');
 
         %mxStep = maxstep(eqdfe_lms,rx_data);
         [eq_y err weights] = eqdfe_lms((rx_data),real(train_symb).');
         
        %[eq_y err weights] = dfe_eq_func(ff,fm0samp,(rx_data),real(train_symb),training_pkt);

        %plot(eq_y(length(train_symb):length(train_symb)+100));

        release(eqdfe_lms);
        pkt_len = 48;

        tr_indx = 1:num_pkt*pkt_len;
        te_indx = tr_indx(end)+1:tr_indx(end)+num_testing_pkts*pkt_len;

        tr_ind = (chunk-off)*pkt_size+1:(chunk-off)*pkt_size+num_pkt*pkt_size;
        te_ind = tr_ind(end)+1:tr_ind(end)+pkt_size*num_testing_pkts;
        %figure;
        bit_vec_tx = reshape(complete_bits(chunk+num_pkt:chunk+num_pkt+ ...
            num_testing_pkts-1,:).',1,prod(size(complete_bits(chunk+ ...
            num_pkt:chunk+num_pkt+num_testing_pkts-1,:))));

        act_data = generate_fm0_sig2(bit_vec_tx,fm0samp).*1e-3;

        %% BER/SNR for test equalizer data
        tot_bits_vec = [];
        bits_dat = eq_y(te_indx);
        %rx_data_dec = rx_data(te_ind);
        rx_data_dec = y_rx(te_ind);
        for y_bit = 1:length(eq_y(te_indx))
           if (bits_dat(y_bit)) < 0
             tot_bits_vec = [tot_bits_vec -0.5.*ones(1,fm0samp/2)]; 
          else
             tot_bits_vec = [tot_bits_vec 0.5.*ones(1,fm0samp/2)];
          end
        end


         tot_bits_a = [];
         tot_bits_b = [];
         act_bits_a = [];
         act_bits_b = [];
         snr_vec =[];
         cnt = 1;
         for kk = 1:num_testing_pkts
              %% Before Equalization
              act_data3 = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
              act_data3 = act_data3(pre_bits+1:end-1);
              pkt = rx_data_dec((kk-1)*pkt_size+1:(kk-1)*pkt_size + pkt_size);
             [snr_basis, dec_bits, sig_power , noise_power, h2, n3] = ...
                 fm0_decode_new_R12(pkt(pre_bits*fm0samp+1-fm0samp:end), ...
                 fm0samp,act_data3);
             tot_bits_b = [tot_bits_b dec_bits];
             bits_b = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
             act_bits_b = [act_bits_b bits_b(pre_bits+1:end-1)];
             ber_b(cnt) = sum(abs(bits_b(pre_bits+1:end-1) - dec_bits))/length(dec_bits);
             snr_vec = [snr_vec snr_basis];
             disp(['Packet BER (Before): ' num2str(ber_b(cnt))]);
             %disp(['Packet SNR in dB   : ' num2str(snr_basis)]);

             %% After Equalization
             pkt = tot_bits_vec((kk-1)*pkt_size+1:(kk-1)*pkt_size + pkt_size).';
             [snr_basis, dec_bits, sig_power , noise_power, h2, n3] = ...
                 fm0_decode_new_R12(pkt(pre_bits*fm0samp+1-fm0samp:end), ...
                 fm0samp,act_data3);
             tot_bits_a = [tot_bits_a dec_bits];
             bits_a = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
             act_bits_a = [act_bits_a bits_a(pre_bits+1:end-1)];
             ber_a(cnt) = sum(abs(bits_a(pre_bits+1:end-1) - dec_bits))/length(dec_bits);
             disp(['Packet BER (After): ' num2str(ber_a(cnt))]);
             cnt =cnt + 1;
         end

        ber_fin_b = sum(abs(tot_bits_b - act_bits_b))/length(tot_bits_b);
        disp(['<strong>Final BER (Before):</strong> ' num2str(ber_fin_b)]);
        ber_fin_a = sum(abs(tot_bits_a - act_bits_a))/length(tot_bits_a);
        disp(['<strong>Final BER (After):</strong> ' num2str(ber_fin_a)]);
        disp(['<strong>Final SNR        :</strong> ' num2str(mean(snr_vec))]);
        ber_vec_a(row,cnttt) = ber_fin_a;
        ber_vec_b(row,cnttt) = ber_fin_b;
        cnttt = cnttt + 1;
        
    end
    row = row + 1;
end


%% weight init

clearvars -except y_rx weights
num_pkt_ind =5;
num_testing_pkts = 300;

W = weights;
%W = 0;
load bit_data
ds_ind = 1;%1;%5;
y_rx = downsample(y_rx,ds_ind);
f_factor = 0.999;%[0.98 0.983 0.989 0.999 1];
%num_pkt_ind = 1;%[90:1:105];
ref_tap = 1;ff_tap = 203;%53;%23; 
fb_tap =100;%200;%100;
const_lvl = max(abs(y_rx(3.79e6:3.86e6))); %(3.126e-4);%3.126e-4; 
%num_pkt_ind = 28*ones(1,20);
Fs = 2e5/ds_ind;
data_rate = 500;%5000;%2e3;
fm0samp = Fs/data_rate;
pkt_bits = 24;
pre_bits = 8;
pkt_size = pkt_bits*fm0samp;


row = 1;
for ff = f_factor
    %clearvars -except ff num_pkt row num_pkt_ind f_factor y_rx complete_bits ber_vec_a ber_vec_b cnttt row 
    
    chunk = 1;
    %num_pkt = 5;
    off = 1;
    temp = chunk;
    training_pkt = [];
    training_bits = [];
    
    cnttt = 1;
    
    for num_pkt = num_pkt_ind
        training_pkt = [];
        training_bits = [];
        temp = chunk;
        for ss = 1:num_pkt
%             if(num_testing_pkts == 26)
%                 disp('ola')
%             end
            training_bits = [training_bits complete_bits(temp,:)];
            training_pkt = [training_pkt generate_fm0_sig2([complete_bits(temp,:)],fm0samp).*const_lvl]; 
            temp = temp + 1;
        end
        
        train_symb = downsample(training_pkt,fm0samp/2);

        rx_data = [y_rx((chunk-off)*pkt_size+1:(chunk-off)*pkt_size+num_pkt*pkt_size+num_testing_pkts*pkt_size)];

         acorr = 0.1;
         eqdfe_lms = comm.DecisionFeedbackEqualizer('Algorithm','RLS','InputSamplesPerSymbol',fm0samp/2, 'ReferenceTap',ref_tap, ...
             'ForgettingFactor',ff,'NumForwardTaps',ff_tap,'AdaptAfterTraining',1,'NumFeedbackTaps',fb_tap,...
             'Constellation',[max(training_pkt) min(training_pkt)],'InitialInverseCorrelationMatrix',acorr,'InitialWeights',(W),'InitialWeightsSource','Property');
 
         %mxStep = maxstep(eqdfe_lms,rx_data);
         [eq_y err weights] = eqdfe_lms((rx_data),real(train_symb).');
         
        %[eq_y err weights] = dfe_eq_func(ff,fm0samp,(rx_data),real(train_symb),training_pkt);

        %plot(eq_y(length(train_symb):length(train_symb)+100));

        release(eqdfe_lms);
        pkt_len = 100;

        tr_indx = 1:num_pkt*pkt_len;
        te_indx = tr_indx(end)+1:tr_indx(end)+num_testing_pkts*pkt_len;

        tr_ind = (chunk-off)*pkt_size+1:(chunk-off)*pkt_size+num_pkt*pkt_size;
        te_ind = tr_ind(end)+1:tr_ind(end)+pkt_size*num_testing_pkts;
        %figure;
        bit_vec_tx = reshape(complete_bits(chunk+num_pkt:chunk+num_pkt+ ...
            num_testing_pkts-1,:).',1,prod(size(complete_bits(chunk+ ...
            num_pkt:chunk+num_pkt+num_testing_pkts-1,:))));

        act_data = generate_fm0_sig2(bit_vec_tx,fm0samp).*1e-3;

        %% BER/SNR for test equalizer data
        tot_bits_vec = [];
        bits_dat = eq_y(te_indx);
        %rx_data_dec = rx_data(te_ind);
        rx_data_dec = y_rx(te_ind);
        for y_bit = 1:length(eq_y(te_indx))
           if (bits_dat(y_bit)) < 0
             tot_bits_vec = [tot_bits_vec -0.5.*ones(1,fm0samp/2)]; 
          else
             tot_bits_vec = [tot_bits_vec 0.5.*ones(1,fm0samp/2)];
          end
        end


         tot_bits_a = [];
         tot_bits_b = [];
         act_bits_a = [];
         act_bits_b = [];
         snr_vec = [];
         cnt = 1;
         for kk = 1:num_testing_pkts
              %% Before Equalization
              act_data3 = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
              act_data3 = act_data3(pre_bits+1:end-1);
              pkt = rx_data_dec((kk-1)*pkt_size+1:(kk-1)*pkt_size + pkt_size);
             [snr_basis, dec_bits, sig_power , noise_power, h2, n3] = ...
                 fm0_decode_new_R12(pkt(pre_bits*fm0samp+1-fm0samp:end), ...
                 fm0samp,act_data3);
             tot_bits_b = [tot_bits_b dec_bits];
             bits_b = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
             act_bits_b = [act_bits_b bits_b(pre_bits+1:end-1)];
             ber_b(cnt) = sum(abs(bits_b(pre_bits+1:end-1) - dec_bits))/length(dec_bits);
             snr_vec = [snr_vec snr_basis];
             disp(['Packet BER (Before): ' num2str(ber_b(cnt))]);
             %disp(['Packet SNR in dB   : ' num2str(snr_basis)]);
             %disp(['Packet BER (Before): ' num2str(ber_b(cnt))]);
             %% After Equalization
             pkt = tot_bits_vec((kk-1)*pkt_size+1:(kk-1)*pkt_size + pkt_size).';
             [snr_basis, dec_bits, sig_power , noise_power, h2, n3] = ...
                 fm0_decode_new_R12(pkt(pre_bits*fm0samp+1-fm0samp:end), ...
                 fm0samp,act_data3);
             tot_bits_a = [tot_bits_a dec_bits];
             bits_a = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
             act_bits_a = [act_bits_a bits_a(pre_bits+1:end-1)];
             ber_a(cnt) = sum(abs(bits_a(pre_bits+1:end-1) - dec_bits))/length(dec_bits);
             disp(['Packet BER (After): ' num2str(ber_a(cnt))]);
             cnt =cnt + 1;
         end

        ber_fin_b = sum(abs(tot_bits_b - act_bits_b))/length(tot_bits_b);
        disp(['<strong>Final BER (Before):</strong> ' num2str(ber_fin_b)]);
        ber_fin_a = sum(abs(tot_bits_a - act_bits_a))/length(tot_bits_a);
        disp(['<strong>Final BER (After):</strong> ' num2str(ber_fin_a)]);
        disp(['<strong>Final SNR        :</strong> ' num2str(mean(snr_vec))]);
        ber_vec_a(row,cnttt) = ber_fin_a;
        ber_vec_b(row,cnttt) = ber_fin_b;
        cnttt = cnttt + 1;
        
    end
    row = row + 1;
end



%num_pkt_ind = [100:100:500, 1000];
%plot(num_pkt_ind,ber_vec_a(1,:),'linewidth',2)
% hold on;plot(num_pkt_ind,ber_vec_a(2,:),'linewidth',2)
%  hold on;plot(num_pkt_ind,ber_vec_a(3,:),'linewidth',2)
%  hold on;plot(num_pkt_ind,ber_vec_a(4,:),'linewidth',2)
%  hold on;plot(num_pkt_ind,ber_vec_a(5,:),'linewidth',2)
% grid on
 %hold on;plot(num_pkt_ind,ber_vec_b(1,:),'linewidth',2)
%legend('0.98', '0.983', '0.989', '0.999', '1','Before Eq');
 %legend('0.999','Before Eq');
% xlabel('Number of testing packets'); ylabel('BER'); 
% title('Distance 4 meters')
