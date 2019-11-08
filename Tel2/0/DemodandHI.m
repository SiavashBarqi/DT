clear
close all
clc

%% Parameters
Ts = 1/1e3;       % symbol time
EbNo = 100000;      % Enrgy of bit over N0
cfo = 1*0.1;      % carrier frequency offset
sfo = .1*1e-3;     % sampling frequency offset
phi = 1*pi/4;    % phase offset
time_offset = 1;  % ideal iz 0=sps=2*sps=...

%% Pulse shaping
a = 0.25;  % roll-off factor
span = 6; % span for pulse shaping
sps = 4;  % sapmple per second
pulseShaping = 1; % 1: RRC 2: Rectang
if pulseShaping == 1
    h = rcosdesign(a,span,sps); % it is normalized energy
    delay = span;
    delay2 = 0; 
    figure(1)
    stem(h,'o')
    legend('Pulse shapping')
elseif pulseShaping == 2
    h = ones(1,sps)/sqrt(sps); % normalized energy
    figure(1)
    stem(h,'o')
    legend('Pulse shapping')
    delay = 1;
    delay2 = sps-1;
end
%% Modulation
M = 8;  % Modulation Order
N = 1e6; % number of symbols
data = randi([0,M-1],1,N);
x = qammod(data, M);
% x = [1 zeros(1,sps*10)]; % for test
%x = pskmod(data, M);
%scatterplot(x)

%% Transmitted Signal
x_upsam = upsample(x,sps);
tx = filter(h,1,x_upsam);
figure (2)
hold on
plot(real(tx(1:sps*min(N,10))))
legend('Real part of transmitted sequence: 10 symbols')

%% Channel and hardware impairements
k = log2(M);
snr = EbNo + 10*log10(k) - 10*log10(sps)
rx_awgn = awgn(tx,snr,'measured');
figure(2)
plot(real(rx_awgn(1:sps*min(N,10))))
legend('Nosiy signal')
hold off

fs = sps/Ts;
figure(3)
hold on
pwelch(rx_awgn,[],[],[],fs,'centered','power')
rx_time_offset = rx_awgn(1+time_offset:end);
t = (1:length(rx_time_offset))/fs;
rx_cfo = rx_time_offset.*exp(1i*2*pi*cfo.*t+1i*phi);
p = 4e4;
q = round(p*(1+sfo));
rx_sfo = resample(rx_cfo,p,q);
rx_HI = rx_sfo;
l = min(min(length(tx),length(rx_awgn)),length(rx_HI));
sig_matrix = [tx(1:l);rx_awgn(1:l);rx_HI(1:l)].';
pwelch(sig_matrix,[],[],[],fs,'centered','power')
legend('noisy recieved signal' , 'data after filtering for transmission'...
    , 'noisy recieved signal' , 'resampled data via inaterpolation')
hold off


%% Matched filter and downsampling at Receiver
rx = filter(h,1,rx_HI);
%stem(rx)
if pulseShaping == 1
    rx_dwnsamp = downsample(rx,sps);
elseif pulseShaping == 2
    rx_dwnsamp = downsample([0 rx],sps); % add zero for rectancular pulses, for RRC we do no need to insert z since the length of filter is sps&span+1
end

figure(4)
ll = 10;
hold on
stem(real(x_upsam(1:ll*sps)),'k')
plot(real(tx(1:ll*sps)))
plot(real(rx(1:ll*sps)./max(real(rx(1:ll*sps)))),'r')
stem(real(upsample(rx(sps:sps:sps*ll),sps)/(max(real(rx(1:ll*sps))))),'*')
legend('tranmsitted bits' , 'Real tranmsitted signal' ,...
    'Real output of matched filter' , 'Real Sampled output')
hold off

figure(5)
hold on
stem(imag(x_upsam(1:ll*sps)),'k')
plot(imag(tx(1:ll*sps)))
plot(imag(rx(1:ll*sps)./max(imag(rx(1:ll*sps)))),'r')
stem(imag(upsample(rx(sps:sps:ll*sps),sps)/(max(imag(rx(1:ll*sps))))),'*')
legend('tranmsitted bits' , 'Imag tranmsitted signal' , ...
    'Imag output of matched filter' , 'Imag Sampled output')
hold off

%%
rx_delay = rx_dwnsamp(delay+1:end);
l = min(length(rx_delay),length(x));
aaa = scatterplot(x(1:l),[],[],'*r');
hold on
scatterplot(rx_delay(1:l),[],[],'.',aaa)
legend('recieved data' , 'data')
hold off

%% BER calculation
Recovered_symnbol = qamdemod(rx_delay,M);
number_of_correct = find(Recovered_symnbol(1:l)-data(1:l)==0);
SER = (l-length(number_of_correct))/l

%% Muller and Muller Algorithm
