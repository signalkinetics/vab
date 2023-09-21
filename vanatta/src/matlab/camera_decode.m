function [bits] = camera_decode(signal,fm0_samps,bit_num)
    s0_1 = transpose(generate_fm0_sig2([1 0 1],fm0_samps));
    s0_1 = s0_1(ceil(fm0_samps/2) + 1 :end - ceil(fm0_samps/2));
    s0_2 = -1.*s0_1;
    
    
    
    s1_1 = transpose(generate_fm0_sig2([1 1 1],fm0_samps));
    s1_1 = s1_1(ceil(fm0_samps/2) + 1 :end - ceil(fm0_samps/2));
    s1_2 = -1.*s1_1;
    
    %seq = [1 0 1 1 0 1 1 1 0 1 0 0 1 0 0 0 1 0 1 0];
    %seq = [0 actual_vec 0];
    %sig_tot = 2*transpose(generate_fm0_sig2(seq,fm0_samps));
        
    bits = [];
    count = 0;

    %signal = signal(ceil(fm0_samps/2) +1 :end-ceil(fm0_samps/2));
    %sig_tot =  sig_tot(ceil(fm0_samps/2) + 1 :end-ceil(fm0_samps/2));
    
    while(count<bit_num)
        
        if (count+2)*fm0_samps  > length(signal)
            end_p = length(signal) ;
        else
            end_p = (count+2)*fm0_samps ;
        end
        
        
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
        if((index==[1]) | (index==[3]))
            bits = [bits, 0];
        end
        if((index==[2]) | (index==[4]))
            bits = [bits, 1];
        end      
        count = count + 1;
    end 
     
    
end