function [y] = my_hpf(x,pass,stop,Fs)

Fpass = pass; %220;%1200;%80;
Fstop = stop; %20;%1000;%20;

xx = x;

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
     %FIRHPF = dsp.HighpassFilter('SampleRate',Fs,...
      %                        'FilterType',filtertype,...
       %                       'PassbandFrequency',Fpass,...
        %                      'StopbandFrequency',Fstop);
    %pad = zeros(chunk*length(xx)/chunk,1);
    pad = zeros(length(xx)/chunk,1);
    temp_hpf = FIRHPF([pad;xx(i:i+length(xx)/chunk-1);pad]);
    temp2_hpf = temp_hpf(length(pad)+del_indx:end-length(pad)+del_indx-1);
    final_hpf = [final_hpf;temp2_hpf];
    
end

y = final_hpf;

end

