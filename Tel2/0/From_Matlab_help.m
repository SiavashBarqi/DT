clear all
clc
%%
EbNo = 7;
rolloff = 0.25; % Filter rolloff
span = 6;       % Filter span
sps = 4;        % Samples per symbol
M = 4;          % Size of the signal constellation
k = log2(M);    % Number of bits per symbol
%%
rrcFilter = rcosdesign(rolloff, span, sps); %discrete sinc
data = randi([0 M-1], 10000, 1);
modData = qammod(data, M);
txSig = upfirdn(modData, rrcFilter, sps); %upsampling
snr = EbNo + 10*log10(k) - 10*log10(sps);
rxSig = txSig + awgn(txSig, snr, 'measured');
rxFilt = upfirdn(rxSig, rrcFilter, 1, sps); %downsampling
rxFilt = rxFilt(span+1:end-span);

h = scatterplot(sqrt(sps)* ... 
    rxSig(1:sps*5000),...
    sps,0,'g.');
hold on;
scatterplot(rxFilt(1:5000),1,0,'kx',h);
%scatterplot(0,1,0,'kx',h);
title('Received Signal, Before and After Filtering');
legend('Before Filtering','After Filtering');
hold off;