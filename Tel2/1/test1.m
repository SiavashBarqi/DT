clear 
clc
close all
%% parameters
snr = 1; % SNR in dB
freq_vec = [5 7 9 13]; % messege frequency vector
fc = min(50, 3*max(freq_vec)); % carrier frequency
ac = 1; %    carrier amplitude
mu = 0.8; % AM modulation index
fs = 3*(fc+max(freq_vec)); % sampling time
t = 0:1/fs:200*1/min(freq_vec); % time vector
%% message gener ation
m_sig = sin(2*pi*freq_vec(1).*t)+2*sin(2*pi*freq_vec(2).*t)+0.8*cos(2*pi*freq_vec(3).*t);
figure(1)
plot(t,m_sig)
hold on
plot(t,-m_sig)
%% carrier
c_sig = ac*cos(2*pi*fc.*t); 
%% modulated signal
% DSB-AM-SC
xc_am_sc = m_sig.*c_sig;
hold on 
plot(t,xc_am_sc,'r')
% DSB-AM
m_tilde = (1+mu*m_sig/max(abs(m_sig)));
xc_am = m_tilde.*c_sig;
figure(2)
hold on
plot(t,m_tilde)
plot(t,xc_am,'r')
% Spectrum analyzer
figure(3)
hold on
pwelch(m_sig,[],[],2^12,fs,'centered','power'); %pwelch(X,WINDOW,NOVERLAP,W,'twosided')
pwelch(m_tilde,[],[],2^12,fs,'centered','power');
pwelch(xc_am_sc,[],[],2^12,fs,'centered','power');
pwelch(xc_am,[],[],2^12,fs,'centered','power');
%% AWGN channel
y_am_sc = awgn(xc_am_sc,snr,'measured');
y_am = awgn(xc_am,snr,'measured');
figure(4)
hold on
pwelch(y_am_sc,[],[],2^12,fs,'centered','power'); %pwelch(X,WINDOW,NOVERLAP,W,'twosided')
pwelch(y_am,[],[],2^12,fs,'centered','power');
%% filtering
y_am_sc_mult = y_am_sc.*c_sig;
pwelch(y_am_sc_mult,[],[],2^12,fs,'centered','power');
Num_filter = [0.000146084040662208,0.000282587801786863,-0.00141681724322881,-0.00319841281840694,0.00622569207782452,0.0176125028881038,-0.0163218970429339,-0.0677229129203946,0.0282733266005913,0.303026235322855,0.466187222473426,0.303026235322855,0.0282733266005913,-0.0677229129203946,-0.0163218970429339,0.0176125028881038,0.00622569207782452,-0.00319841281840694,-0.00141681724322881,0.000282587801786863,0.000146084040662208];%[0.00123364867740107,0.000742224511005134,-0.00588161637091575,-0.00594255931827098,0.0157865984257343,0.0242729804262498,-0.0297401543465135,-0.0760474411726777,0.0425486702035788,0.306958406448995,0.452207574827808,0.306958406448995,0.0425486702035788,-0.0760474411726777,-0.0297401543465135,0.0242729804262498,0.0157865984257343,-0.00594255931827098,-0.00588161637091575,0.000742224511005134,0.00123364867740107];
am_sc_out = filter(Num_filter,1,y_am_sc_mult);
figure(5)
plot(am_sc_out)
hold on
plot(m_sig)
pause(3)
close all
% Spectrum
figure(6)
hold on
pwelch(am_sc_out,[],[],2^12,fs,'centered','power');
pwelch(m_sig,[],[],2^12,fs,'centered','power');
%%
sound(m_sig)
sound(am_sc_out)

%% FM Modulator and Demodulator
kf = 50; % modulation index or frequency deviation
kp = 50; % PM mod. index
m_sig = sin(2*pi*freq_vec(1).*t);%+sin(2*pi*freq_vec(2).*t)+cos(2*pi*freq_vec(3).*t);
m_sig_nrm = m_sig/max(abs(m_sig)); % normalized message
x_pm = cos(2*pi*fc*t+kp*m_sig_nrm);
figure(7)
pwelch(x_pm,[],[],2^12,fs,'centered','power');
