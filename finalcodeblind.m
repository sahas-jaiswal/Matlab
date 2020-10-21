% This code is to plot receiver operating characteristic curve for simple energy
% detection, when the primary signal is real Gaussian signal and noise is
% addive white real Gaussian. Here, the threshold is available
% analytically.
% Code written by: Sanket Kalamkar, Indian Institute of Technology Kanpur,
% India.


clc
close all
clear all
L = 1000;
snr_dB = -10; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf)
    m
    i = 0;
for kk=1:10000 % Number of Monte Carlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 y = s + n; % Received signal at SU
 energy = abs(y).^2; % Energy of received signal over N samples
 energy_fin =(1/L).*sum(energy); % Test Statistic for the energy detection
 thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1; % Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
 if(energy_fin >= thresh(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i = i+1;
 end
end
Pd(m) = i/kk; 
end
plot(Pf, Pd)
hold on
%% Theroretical ecpression of Probability of Detection; refer above reference.
thresh = (qfuncinv(Pf)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
plot(Pf, Pd_the, 'r')
hold on
% Step 1. Generate message signal of length >= M^L.
M = 16; % Alphabet size of modulation
L = 1; % Length of impulse response of channel
msg = [0:M-1 0]; % M-ary message sequence of length > M^L

% Step 2. Modulate the message signal using baseband modulation.
hMod = comm.RectangularQAMModulator(M);  % Use 16-QAM.
modsig = step(hMod,msg'); % Modulate data
Nsamp = 16;
modsig = rectpulse(modsig,Nsamp); % Use rectangular pulse shaping.

% Step 3. Apply a transmit filter.
txsig = modsig; % No filter in this example

% Step 4. Run txsig through a noiseless channel.
rxsig = txsig*exp(1i*pi/180); % Static phase offset of 1 degree
% Step 5. Use the semianalytic function.
% Specify the receive filter as a pair of input arguments.
% In this case, num and den describe an ideal integrator.
num = ones(Nsamp,1)/Nsamp;
den = 1;
EbNo = 0:20; % Range of Eb/No values under study
ber = semianalytic(txsig,rxsig,'qam',M,Nsamp,num,den,EbNo);

% For comparison, calculate theoretical BER.
bertheory = berawgn(EbNo,'qam',M);

% Plot computed BER and theoretical BER.
figure; semilogy(EbNo,ber,'k*');
hold on; semilogy(EbNo,bertheory,'ro');
title('Semianalytic Communication error Compared with old False?alarm prababality');
legend('Semianalytic BER with Phase Offset',...
   ' error Without Phase Offset','Location','SouthWest');
hold off;
