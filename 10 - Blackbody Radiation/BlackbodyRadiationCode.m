%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% Blackbody Radiation                                       %                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment setup:
Ad = 1*10^(-6); % Detector Area, m^2

sourceAngle = 0; % Both the collector and emitter are flat
sensorAngle = 0;

% Inches are converted to meters, C to K, qs are in Watts

% Values when held constant:
D = 0.6*0.0254;     % Diameter, in to m
T = 605.9 + 273.15; % Temperature, C to K
d = 9*0.0254;       % Distance from aperature

% Trial 1: Constant temperature and distance to aperature
D1 = [0.6000, 0.4000, 0.2000, 0.1000, 0.0500]*0.0254;
q1 = [  30.2,   13.9,   3.45,   0.72,   0.13]*10^(-6);
T1 = [     T,      T,      T,      T,      T];
d1 = [     d,      d,      d,      d,      d];

% Trial 2: Constant aperature diameter and temperature
D2 = [   D,    D,    D,    D,    D];
q2 = [3.87, 3.04, 2.31, 1.91, 1.61]*10^(-5);
T2 = [   T,    T,    T,    T,    T];
d2 = [   8,    9,   10,   11,   12]*0.0254;

% Trial 3: Constant distance and aperature diameter 
D3 = [   D,    D,    D,    D,    D,    D,    D,    D,    D,    D];
q3 = [3.00, 3.37, 3.67, 3.97, 4.24, 4.48, 4.70, 4.95, 5.26, 5.55]*10^(-5);
T3 = [ 606,  620,  635,  650,  665,  680,  695,  710,  725,  740] + 273.15;
d3 = [   d,    d,    d,    d,    d,    d,    d,    d,    d,    d];

% All trials
DAll = [D1, D2, D3];
qAll = [q1, q2, q3];
TAll = [T1, T2, T3];
dAll = [d1, d2, d3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate radiative heat transfer based on experiment parameters
Eb = (5.67e-8)*TAll.^4;   % Blackbody emissive power
Ib = Eb/pi;               % Diffuse Blackbody emission intensity
As = pi*(DAll.^2)/4;      % Area of source
An = Ad*cos(sensorAngle); % Normal detector area
solidAngle = An./dAll.^2; % Solid angle of detector

qCalc = Ib.*As.*cos(sourceAngle).*solidAngle; % Radiative heat transfer

% Plot measured vs calculated heat transfer, overlay a 1:1 line, annotate
figure; hold on
plot(qCalc, qAll, 'kx');
xl = xlim;
OneToOne = linspace(xl(1), xl(2), 100);
plot(OneToOne, OneToOne, 'b');
xlabel('Measured radiative heat transfer (Watts)');
ylabel('Calculated radiative heat transfer (Watts)');
legend('Data', '1:1 Line', 'location', 'southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate blackbody temperature
Ib2 = qAll ./ (As.*cos(sourceAngle).*solidAngle); % Blackbody intensity
Eb2 = pi*Ib2;                                     % Emissive power
TCalc = (Eb2 / (5.67e-8)).^(1/4);                 % Blackbody temperature

% Plot measured vs calculated temperature, overlay a 1:1 line, annotate
figure; hold on
plot(TCalc, TAll, 'kx');
xl = xlim;
OneToOne = linspace(xl(1), xl(2), 100);
plot(OneToOne, OneToOne, 'b');
xlabel('Measured body temperature (K)');
ylabel('Calculated body temperature (K)');
legend('Data', '1:1 Line', 'location', 'southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold on
plot(T3.^4, q3, 'kx');

lfit3 = polyfit(T3.^4, q3, 1);
plot(T3.^4, lfit3(1)*T3.^4 + lfit3(2), 'k-');

xlabel('Temperature^4 T^4 (K^4)');
ylabel('Measured radiative heat transfer (W)');
fiteq = sprintf('q = %.2s*T^4 + %.2s', lfit3(1), lfit3(2));
legend('Data', fiteq, 'location', 'southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold on
plot(1./d2.^2, q2, 'kx');

lfit4 = polyfit(1./d2.^2, q2, 1);
plot(1./d2.^2, lfit4(1)./d2.^2 + lfit4(2), 'k-');

xlabel('Inverse of distance squared 1/r^2 (m^-2)');
ylabel('Measured radiative heat transfer (W)');
fiteq = sprintf('q = %.2s*1/r^2 + %.2s', lfit4(1), lfit4(2));
legend('Data', fiteq, 'location', 'southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold on
plot(D1.^2, q1, 'kx');

lfit5 = polyfit(D1.^2, q1, 1);
plot(D1.^2, lfit5(1)*D1.^2 + lfit5(2), 'k-');

xlabel('Aperature Diameter Squared D^2 (m^2)');
ylabel('Measured radiative heat transfer (W)');
fiteq = sprintf('q = %.2s*D^2 + %.2s', lfit5(1), lfit5(2));
legend('Data', fiteq, 'location', 'southeast');

