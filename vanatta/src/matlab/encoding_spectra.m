fs = 2e5;
N_bits = 1024;

fb = 250;
Tb = 1/fb;

code = randi([0 1],1,N_bits); 

%% RECT PULSE (NRZ CODE)
nrz_seq = repelem(code,fs/fb);
seq = nrz_seq;

%% FM0
fm0_seq = generate_fm0_sig2(code,fs/fb);
t = [0:1/fs:N_bits/fb-1/fs];

seq = fm0_seq;

%% MILLER
m=2;
blf = fb*m;
miller_seq = miller_encode(code,m,blf,fs);

lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',blf/4 ...
                    ,'StopbandFrequency',blf/2,'StopbandAttenuation' ...
                    ,30,'PassbandRipple',1,'DesignMethod','kaiserwin','SampleRate',fs);

% bpFilt = designfilt('bandpassfir' ...
%                     ,'PassbandFrequency1',m/2*blf-2*blf/m ,'StopbandFrequency1',(0.5)*(m/2*blf-2*blf/m) ... 
%                    ,'StopbandAttenuation1',30 ...
%                    ,'PassbandFrequency2',m/2*blf+2*blf/m ,'StopbandFrequency2',(1.5)*(m/2*blf+2*blf/m) ... 
%                    ,'StopbandAttenuation2',30, ...
%                    'PassbandRipple',1,'DesignMethod','kaiserwin','SampleRate',fs);

t = [0:1/fs:(length(miller_seq)-1)/fs];

seq = miller_seq;

lo = sin(2*pi*blf*t);

seq = miller_seq.*lo;

seq = filtfilt(lpFilt,seq);


%% PLOTTING


t_window = length(t)/10/fs;
window = chebwin(t_window*fs);

[pxx,f] = pwelch(seq,window,[],[],fs);

[pxx_orig,f] = pwelch(miller_seq,window,[],[],fs);

figure(1);
subplot(3,1,1);
%hold on
plot(t,seq);
hold on;
% plot(t,miller_seq);
grid on;
grid minor;
xlabel("Time (s)");
ylabel("Value");
ylim([-1.5 1.5]);
subplot(3,1,2);
hold on
plot(f/1e3,10*log10(pxx));
% ylim([-50 0]);
xlim([0 10]);
xlabel("Freq (kHz)");
ylabel("Magnitude (dB20)");
grid on;
grid minor;

subplot(3,1,3);

% plot(f/1e3,10*log10(pxx_orig));
% % ylim([-50 0]);
% xlim([0 10]);
% xlabel("Freq (kHz)");
% ylabel("Magnitude (dB20)");
% grid on;
% grid minor;

function out = rect(t)
    out = zeros(size(t),'like',t);

    out(logical((t < 1/2) .* (t > -1/2))) = 1;
end

function out = miller_encode(code,m,blf,fs)
    s1 = [1 1];
    s2 = [1 -1];
    s3 = -s2;
    s4 = -s1;

    if rem(log2(m),1) ~= 0
        error("M must be a power of 2");
    end

    data_mat = [s1;s2;s3;s4];

    T = [4 2;4 3;1 2;1 3];

    out = zeros(1,2*length(code));

    if code(1) == 0
        out(1:2) = s1;
        state = 1;
    elseif code(1) == 1
        out(1:2) = s2;
        state = 2;
    end

    for i=2:length(code)
        bit = code(i);
        state = T(state,bit+1);
        out(2*i-1:2*i) = data_mat(state,:);
    end
    
    if m ~= 0
        out = repelem(out,m).*repmat([1 -1],1,floor(length(out)*m/2));
    end

    out = repelem(out,fs/blf/2);
end