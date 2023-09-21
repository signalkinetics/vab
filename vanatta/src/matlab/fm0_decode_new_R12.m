function [snr_basis, bits, sig_power , noise_power, chan_vec, n_vec] = fm0_decode_new_R12(signal,fm0_samps,actual_vec)
    s0_1 = transpose(generate_fm0_sig2([1 0 1],fm0_samps));
    s0_1 = s0_1(ceil(fm0_samps/2) + 1 :end - ceil(fm0_samps/2));
    s0_2 = -1.*s0_1;
    
    
    s1_1 = transpose(generate_fm0_sig2([1 1 1],fm0_samps));
    s1_1 = s1_1(ceil(fm0_samps/2) + 1 :end - ceil(fm0_samps/2));
    s1_2 = -1.*s1_1;
    

    
    %seq = [1 0 1 1 0 1 1 1 0 1 0 0 1 0 0 0 1 0 1 0];
    seq = [0 actual_vec 0];
    sig_tot = 2*transpose(generate_fm0_sig2(seq,fm0_samps));
        
    bits = [];
    count = 0;

    signal = signal(ceil(fm0_samps/2) +1 :end-ceil(fm0_samps/2)+1);
    sig_tot2 =  sig_tot(ceil(fm0_samps/2) + 1 :end-ceil(fm0_samps/2)+1);
    sig_tot = sig_tot2;
    while(count<length(actual_vec))
        
        if (count+2)*fm0_samps  > length(signal)
            end_p = length(signal) ;
        else
            end_p = (count+2)*fm0_samps ;
        end
        count;
        
        new_signal = signal(count*fm0_samps + 1 :end_p);
        
        
        new_signal = new_signal - mean(new_signal);
        s0_1 = s0_1 - mean(s0_1);
        s1_1 = s1_1 - mean(s1_1);
        
        new_signal = new_signal(ceil(fm0_samps/2) + 1:end-ceil(fm0_samps/2));
        s0_1t =  s0_1(ceil(fm0_samps/2) + 1:end-ceil(fm0_samps/2));
        s1_1t =  s1_1(ceil(fm0_samps/2) + 1:end-ceil(fm0_samps/2));
        s0_2t = s0_2(ceil(fm0_samps/2) + 1:end-ceil(fm0_samps/2));
        s1_2t = s1_2(ceil(fm0_samps/2) + 1:end-ceil(fm0_samps/2)); 
        corr01 = (1/norm(s0_1t).^2)*sum((new_signal).' * (s0_1t));
        corr11 = (1/norm(s1_1t).^2)*sum((new_signal).' * (s1_1t));
        corr02 = (1/norm(s0_2t).^2)*sum((new_signal).*(s0_2t));
        corr12 = (1/norm(s1_2t).^2)*sum((new_signal).*(s1_2t));

        four_corr = [abs(corr01), abs(corr11),abs(corr02), abs(corr12)];
        
        %proj(count+1,:) = four_corr;
        
        index = find(four_corr==max(four_corr));
        index;
        if((index==[1]) | (index==[3]))
            bits = [bits, 0];
        
        else if((index==[2]) | (index==[4]))
            bits = [bits, 1];
        else 
            bits=[bits, 0];
        end
        end
        count = count + 1;
    end 
    length(actual_vec);
    %% SNR calculation
    
    %% PART1 : packet h calculation:
    cnt = 1;
    count = 0;
    h_tot = [];
    pol= [];
    while (count < length(actual_vec))
        
        if (count+2)*fm0_samps > length(signal)
            end_p = length(signal);
        else
            end_p = (count+2)*fm0_samps;
        end
        
        
        new_signal = signal(count*fm0_samps + 1 :end_p);

        act_signal = sig_tot(count*fm0_samps+1:end_p);
        
        new_signal = new_signal - mean(new_signal);
        act_signal = act_signal - mean(act_signal);
        
        
        new_signal = new_signal(ceil(fm0_samps/2)+1:end-ceil(fm0_samps/2));
        act_signal = act_signal(ceil(fm0_samps/2)+1:end-ceil(fm0_samps/2));
        
        
    
        h_tot(cnt) = (1/norm(act_signal).^2).*sum(new_signal.*act_signal);

        
        count = count + 1;
        cnt = cnt+1;
    end
    
    %h = 1 ; %mean(h_tot); %% Hardcoded Channel
    %h_tot=mean(h_tot)*ones(20,1);
%     pol(1,:) = ones(1,length(pol)); %% Hardcoded Polarity
    %% PART2: packer snr calculation
    
    cnt = 1;
    count = 0;
    snr2 = [];
    n_tot = [];
    while (count < length(actual_vec))
        
        if (count+2)*fm0_samps > length(signal)
            end_p = length(signal);
        else
            end_p = (count+2)*fm0_samps;
        end
        
        
        new_signal = signal(count*fm0_samps+1:end_p);

        act_signal = sig_tot(count*fm0_samps+1:end_p);
        
        new_signal = new_signal - mean(new_signal);
        act_signal = act_signal - mean(act_signal);
        
        
        new_signal = new_signal(ceil(fm0_samps/2)+1:end-ceil(fm0_samps/2));
        act_signal = act_signal(ceil(fm0_samps/2)+1:end-ceil(fm0_samps/2));
        
        
    
        %h(cnt) = (1/length(act_signal)).*sum(new_signal.*act_signal);
        %h(cnt) = (1/norm(act_signal).^2).*sum(new_signal.*act_signal);
        
         bits_2 = seq(count+1);
         bits_inv = ~bits_2;
         n_basis = 2*transpose(generate_fm0_sig2([bits_inv],fm0_samps));
%         
         bits_2 = seq(count+1);
         bits_inv = bits_2;
         n_basis2 = 2*transpose(generate_fm0_sig2([bits_inv],fm0_samps));
        
        
        %n_basis = n_basis(ceil(fm0_samps/2)+1:end - ceil(fm0_samps/2));
        
        %n1 = 1/sqrt(sum(n_basis(1:end/2).*n_basis(1:end/2)))*sum(new_signal(1:end/2).*n_basis(1:end/2));
        %n2 = 1/sqrt(sum(n_basis(end/2+1:end).*n_basis(end/2+1:end)))*sum(new_signal(end/2+1:end).*n_basis(end/2+1:end));
        
        %n =  [1/sqrt(sum(n_basis(1:end/2).*n_basis(1:end/2))).*n1.*n_basis(1:end/2); zeros(size(n_basis(end/2+1:end)))] ...
        %      +  [zeros(size(n_basis(1:end/2))); 1/sqrt(sum(n_basis(end/2+1:end).*n_basis(end/2+1:end))).*n2.*n_basis(end/2+1:end)];
       
         %n = new_signal- pol(cnt)*sqrt(h).*act_signal;
         n = new_signal -  h_tot(cnt).*act_signal;
         n_tot = [n_tot; n];

         
         %%%%%% Saad Added This %%%%%%%%%%
%          
           nb = n(1:end);
           nb = (1/(norm(n_basis).^2)).*sum(n_basis.*nb).*n_basis;
           % (1/(norm(n_basis2).^2)).*sum(n_basis2.*nb).*n_basis2;
         
         %n = mean(n(1:ceil(fm0_samps/2)).^2) + mean(n(ceil(fm0_samps/2)+1:end).^2);
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%	%%
         
         %n = n(1:end-1);
         %n = (1/(norm(n_basis).^2)).*sum(n_basis.*n).*n_basis;
         
         %n = mean(n(1:ceil(fm0_samps/2)).^2) + mean(n(ceil(fm0_samps/2)+1:end).^2);
         


%         n_half1 = (1/norm(n_basis(1:end/2)).^2).*sum(n_basis(1:end/2).*n(1:end/2)).*n_basis(1:end/2);
%         n_half2 = (1/norm(n_basis(end/2+1:end)).^2).*sum(n_basis(end/2+1:end).*n(end/2+1:end)).*n_basis(end/2+1:end);
%         n = [n_half1; n_half2];
        
        
        
        %nn(cnt) = (1/norm(act_signal).^2).*sum(abs(n).^2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %nn(cnt) = (1/length(n).*sum(n.^2));
        %nn(cnt) = mean(n.^2); % normalize the noise to get SNR properly computed

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         h(cnt)=sum(new_signal.*act_signal).*((1/sqrt(sum(act_signal.*act_signal))));
%         n = new_signal - h(cnt).*act_signal;
%         nn(cnt) = (sum(abs(n).^2));
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        %snr(cnt) = (abs(h(cnt))^2)/nn(cnt);
        %snr2(cnt) = snr(h(cnt).*act_signal(1:end-1)+n,n);
        
        %h = h/sqrt(2);
        
        %snr2(cnt) = 10*log10((sum((h.*(act_signal.^2))))/((sum(n.^2)))); % originally : mean
        
        %snr2(cnt) = 10*log10((sum((abs(h*(act_signal)))))/((sum(n.^2)))); % originally : mean
%         snr2(cnt) = 10*log10((sum((abs((h_tot(cnt)*act_signal).^2))))/((sum(n.^2)))); % originally :
        
%         snrb(cnt) =  10*log10((sum((abs((h_tot(cnt)*act_signal).^2))))/((sum(nb.^2))));
        
        
        %snr2(cnt) = 10*log10((sum((h(cnt)).^2))/sum(n.^2));
        %snr3=10*log10(h(cnt)^2/mean(n.^2));
        %snr2(cnt) = (h(cnt)./sqrt(2)).^2/n;
        h2(cnt) = mean((abs(h_tot(cnt)*(act_signal)).^2));
        
        n2(cnt) = mean(abs(n).^2);
        n3(cnt) = mean(abs(nb).^2);
        %snr5(cnt) = h(cnt).^2/var(n);
        %snr3(cnt) =  sum((h(cnt).*act_signal).^2)/nn(cnt);
        count = count + 1;
        cnt = cnt+1;
    end
    
%      snr_mean =  mean(snr2) ; %10*log10(mean(snr2)); ;%
     
     hhh = mean(abs(h2));
     chan_vec = h2;
     n_vec = n3;
     nnn = mean(n2);
     nnn2 = mean(n3); % noise using basis 
%      snrr=mean(hhh/nnn);
     snr_mean = 10*log10(hhh/nnn);% + 10*log10(fm0_samps);
     snr_basis = abs((hhh-nnn2)/(2*nnn2)); %10*log10(hhh/nnn2);
     sig_power = abs(hhh-nnn2);
     noise_power = 2*nnn2;
%      
%      snr_basis = mean(snrb);
%     
end