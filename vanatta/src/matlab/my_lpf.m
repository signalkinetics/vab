function [y] = my_lpf(x,pass,stop,Fs)

% Low pass
FIRLPF = dsp.LowpassFilter('SampleRate',Fs,'PassbandFrequency',pass, 'StopbandFrequency',stop);

shifted_sig = FIRLPF([1 zeros(1,Fs)].');
[max_filt,del_indxlp] = max(shifted_sig);
disp(['Low Pass Filter Delay(Samples): ' num2str(del_indxlp)]);
release(FIRLPF)

    final_lpf = [];
    chunk = 1;
    mult_s = x.';
    for i = 1:length(mult_s)/chunk:length(mult_s)

       pad = zeros(1*length(mult_s)/chunk,1);
       temp_lpf = FIRLPF([pad;mult_s(i:i+length(mult_s)/chunk-1);pad]);
       temp2_lpf = temp_lpf(length(pad)+del_indxlp:end-length(pad)+del_indxlp-1);
       final_lpf = [final_lpf;temp2_lpf];
       i;

    end
    
y = final_lpf;

end

