clear; %Clear all variables in the workspace 
close all; %Closes all figure windows
clc %Clear command window

%Load the data for Assignment 1, which 
%contains the variables ComplexAmplitude and TxRxDist:
load ('Assignment1'); 

RxPow=(abs(ComplexAmplitude).^2); %Received power

%(note: the Tx power has been removed, so the data is in dB, not in dBm!)

Dist_dB=10*log10(TxRxDist); %Tx-Rx separation in dB.
 
D_win=[]; P_SSA=[]; ssfamplitude=[]; %preallocation of variables in the for-loop
for i=6:5:length(ComplexAmplitude)-5
     D_win=[D_win; mean(TxRxDist(i-5:i+5))]; %Average TxRxDistance for the small-scale averaged power
     P_SSA=[P_SSA; mean(RxPow(i-5:i+5))]; %Small scale averaged power.
     
     ratio=(abs(ComplexAmplitude(i-5:i+5)))./(mean(abs(ComplexAmplitude(i-5:i+5))));
     ssfamplitude=[ssfamplitude ratio]; %Small-scale fading amplitude based on the small-scale averaged power. 
end


figure(1) %new figure window
semilogx((TxRxDist),10*log10(RxPow),'k.') %plot distance against received power, with logarithmic x-axis.
xlabel('Tx-Rx separation distance [m]', 'FontSize', 13, 'Interpreter', 'latex') %label for the x-axis
ylabel('Channel gain [dB]', 'FontSize', 13, 'Interpreter', 'latex') %label for the y-axis.
grid on %add grid
hold on %hold on, to be able to plot another figure, while keeping the old plot.
semilogx((D_win),10*log10(P_SSA),'bo') %Plot small-scale averaged channel gain

%A plot of the free-space pathloss at 2.6 GHz, using the same figure as above:
semilogx(([10 200]),-20*log10(4*pi*[10 200]/(3e8/2.6e9)),'r-')
%A legend of for window 1:
legend('Channel gain','Small-scale averaged channel gain','Free space channel gain', 'Location', 'best', 'Location', 'best', 'Interpreter','latex', 'Fontsize', 11)

%--------------------------------------------------

%A2: Find the ordinary least-squares estimate of PL(d0)
%and pathloss exponent n:
d_0 = 1;
fun = @(x,D_win)x(1)-10*x(2)*log10(D_win/d_0); %P(d) = P(d_0) - 10*n*log(10/d_0)
x = lsqcurvefit(fun,[0,0],D_win,10*log10(P_SSA));

PL_d = x(1) - 10*x(2)*log10(D_win/d_0);

%Plot the average received power using the pathloss estimates found above: 
figure(2)
xlabel('Tx-Rx separation distance [m]', 'FontSize', 13, 'Interpreter', 'latex') %label for the x-axis
ylabel('Channel gain [dB]', 'FontSize', 13, 'Interpreter', 'latex') %label for the y-axis.
grid on %add grid
hold on %hold on, to be able to plot another figure, while keeping the old plot.
plot((D_win),10*log10(P_SSA),'bo',D_win,fun(x,D_win),'g-', 'LineWidth', 1.2)
legend('Small-scale averaged channel gain','OLS fitted curve', 'Location', 'best', 'Location', 'best', 'Interpreter','latex', 'Fontsize', 11)

%--------------------------------------------------

%A3: In a new window, plot the empirical cdf of the large-scale fading from
%the measurement:
figure(3)

LSF = P_SSA - PL_d;

%Plot the empirical cumulative distribution funcfion
e_LSF = cdfplot(LSF); e_LSF.LineStyle = 'none'; e_LSF.Marker = '.';
hold on; grid on; 

%Plot the CDF for the normal distribution based on the maximum-likelihood estimate
%for a normal distribution: 
N_LSF = length(LSF);

mean_LSF = 1/N_LSF * sum(LSF);
var_LSF = 1/N_LSF * sum((LSF - mean_LSF).^2);

CDF_Normal_LSF = (1/2) * (1+erf((LSF-mean_LSF)/(var_LSF*2^(-2))));
plot(LSF, CDF_Normal_LSF, '.')
legend('Emprical', 'Gaussian Modelled', 'Location', 'best', 'Interpreter','latex', 'Fontsize', 11)
xlabel('Large Scale Fading [dB]', 'FontSize', 13, 'Interpreter', 'latex') %label for the x-axis
ylabel('cumulative distribution function of LSF', 'FontSize', 13, 'Interpreter', 'latex')%, xlim([min(LSF) max(LSF)]) %label for the y-axis.
title('\textbf{CDF of Large Scale Fading}', 'Interpreter', 'latex', 'FontSize', 14)
%--------------------------------------------------

%A4: In a new window, plot the empirical cdf of the small scale fading amplitude
figure(4)
e_SSF = cdfplot(ssfamplitude); e_SSF.LineStyle = 'none'; e_SSF.Marker = '.';
hold on; grid on;
%--------------------------------------------------
%Use the maximum-likelihood expression for the
%Rayleigh distribution to find an estimate of the scale
%parameter sigma:
N_SSF = length(ssfamplitude);

var_SSF = 1/(2*N_SSF) * sum(ssfamplitude.^2); % equation (6)
%--------------------------------------------------
%Plot the cdf of the Rayleigh distribution, using
% the estimate obtained in A4:
CDF_Ray_SSF = 1 - exp(-(ssfamplitude.^2)/(2*var_SSF)); % equation (7)
plot(ssfamplitude, CDF_Ray_SSF, '.')
legend('Emprical', 'Rayleigh Modelled', 'Location', 'best', 'Interpreter','latex', 'Fontsize', 11)
xlabel('Small Scale Fading amplitude', 'FontSize', 13, 'Interpreter', 'latex') %label for the x-axis
ylabel('cumulative distribution function of $SSF_{amp}$', 'FontSize', 13, 'Interpreter', 'latex') %label for the y-axis.
title('\textbf{CDF of Small Scale Fading Amplitude}', 'Interpreter', 'latex', 'FontSize', 14)
%-------------------------------------------------- 


%--------------------------------------------------
%A5: Plot the cdf for the Rice distribution with the
%parameter estimates given in the assignment. Hint:
%Type "help marcumq".
sigRice=0.489;
v=0.84185;

CDF_Rice_SSF = 1 - marcumq(v/sigRice, ssfamplitude/sigRice);

figure(5)
e_SSF = cdfplot(ssfamplitude); e_SSF.LineStyle = 'none'; e_SSF.Marker = '.';
hold on; grid on;
plot(ssfamplitude, CDF_Ray_SSF, '.')
plot(ssfamplitude, CDF_Rice_SSF, '.');
legend('Emprical', 'Rayleigh Modelled', 'Rice Modelled', 'Location', 'best', 'Interpreter', 'latex', 'Fontsize', 11)
xlabel('Small Scale Fading amplitude', 'FontSize', 13, 'Interpreter', 'latex') %label for the x-axis
ylabel('cumulative distribution function of $SSF_{amp}$','FontSize', 13, 'Interpreter', 'latex') %label for the y-axis.
title('\textbf{CDF of Small Scale Fading Amplitude}', 'Interpreter', 'latex', 'FontSize', 14)
%-------------------------
