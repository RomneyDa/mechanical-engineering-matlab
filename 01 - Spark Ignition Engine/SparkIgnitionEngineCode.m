% Dallin Romney
% Spark Ignition Engine

clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%
Speed = [2995, 2794, 2598, 2396, 2197]; % Engine Speed, RPM
Torque = [2.21, 2.74, 2.93, 3.10, 3.01]; % Engine Torque, Nm
Fuel = [2.00, 2.00, 2.00, 2.00, 2.00]; % Cm (5.1 mL/cm)
Time = [32.9, 37.5, 41.2, 42.3, 44.5]; % Time, s
Tfuel = [21.7, 21.9, 22.0, 22.1, 22.2]; % Fuel temperature, deg C
Tinlet = [19.8, 20.0, 20.1, 20.2, 20.2]; % Inlet temperature, deg C
AirCon = [177, 168, 166, 165, 163]; % Air consumption rate, L/min
Texhau = [505, 490, 481, 460, 446]; % Exhaust temperature, deg C

r = 7.1; % Compression ratio
k = 1.33; % Ratio of specific heats
Troom = 21.3; % Room temp, deg C
v = 0.148; % Swept volume in L
P = 867.1; % Room pressure (mBar)
LHV = 44000; % Lower heating value, kJ/kg
dFuel = 726; % Fuel density, kg/m^3
DrySpeed = [2514, 2114, 2805, 1517]; % RPM
DryTorque = [-1.92, -1.88, -1.79, -1.76]; % Nm

%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%

% Engine output power in kW and horsepower
SpeedRad = Speed*pi/30; % Convert Speed to rad/s
Power = Torque.*SpeedRad/1000; % Engine output power, kW
PowerHP = Power*1.34102; % Engine output power, hp

% Rate of energy input to the engine in the fuel in kW and horsepower
% 1 cm drop in flow tube = 5.1 cm^3 or mL
mFuel = Fuel*5.1/10^6*dFuel./Time; % Mass flow rate in of fuel, kg/s
Ein = mFuel*LHV; % Energy in, kW

% Rate of mechanical losses due to friction and inertia in kW and horsepower
% at each engine speed RPM using linear interpolation and or linear
% extrapolation.
DrySpeedRad = DrySpeed*pi/30;
fLossDry = -DrySpeedRad.*DryTorque/1000;
p = polyfit(DrySpeedRad, fLossDry, 1);
fLoss = SpeedRad*p(1) + p(2); % Frictional losses, kW
fLosshp = fLoss*1.34102; % Frictional losses, hp

% Miscellaneous loss rate in kW and HP
mLoss = -(Power + fLoss - Ein); % Miscellaneous Losses, kW
mLosshp = mLoss*1.34102; % Miscellaneous Losses, hp

% Thermal efficiency in %
n = Power./Ein;

% Mean effective pressure in kPa and psi
MEP = 4*pi*Torque/(1000*v); % Mean Effective Pressure [kPa]
MEPpsi = MEP/7; % Mean Effective Pressure [psi]

% Ideal Otto cycle thermal efficiency for the engine assuming the ratio of
% specific heats k = 1.33
nOtto = 1 - 1/(r^(k-1)); % Ideal otto engine efficiency

%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%

% 1.a. Torque [Nm] vs. Engine speed [RPM]
figure; plot(Speed, Torque);
title('Plot 1.a: Torque vs. Speed');
xlabel('Speed [RPM]'); ylabel('Torque [Nm]');

% 1.b. Power output [kW] versus engine speed [RPM] (data & theoretical)
figure; plot(Speed, Power);
title('Plot 1.b: Power vs. Speed');
xlabel('Speed [RPM]'); ylabel('Power Output [kW]');
legend('data', 'theoretical');

% 1.c. Fuel flow rate / useful power [kg*kW/s] versus engine speed [RPM]
figure; plot(Speed, mFuel./Power);
title('Plot 1.c: Useful power vs. Speed');
xlabel('Speed [RPM]'); ylabel('Useful power [kg*kW/s]');

% 1.d. Thermal efficiency [%] versus engine speed in RPM (data + theor)
figure; hold on;
plot(Speed, n*100); plot(Speed, nOtto*100);
title('Plot 1.d: Thermal Efficiency vs. Speed');
xlabel('Speed [RPM]'); ylabel('Thermal Efficiency [%]');
legend('data', 'theoretical');

% 1.e. Rate of fuel energy input (open circles), frictional and inertial
% losses (open squares), and miscellaneous loss rate (open diamonds) [kW]
% vs. engine speed [RPM]
figure; hold on;
plot(Speed, Ein, 'o'); plot(Speed, fLoss, 's'); plot(Speed, mLoss, 'd');
title('Plot 1.e: Fuel Energy Input and Losses vs. Speed');
xlabel('Speed [RPM]'); ylabel('[kW]');
legend('Fuel Energy Input', 'Frictional/Inertial Losses', 'Miscellaneous Loss Rate');

% 1.f. Mean Effective Pressure [kPa] versus engine speed [RPM]
figure; plot(Speed, MEP);
title('Plot 1.f: Mean Effective Pressure vs. Speed');
xlabel('Speed [RPM]'); ylabel('Mean Effective Pressure [kPa]');

%%%%%%%%%%%%%%%% SHORT ANSWER QUESTIONS %%%%%%%%%%%%%%%%
%
% 2a. Discuss the percent contribution of each of the terms in the energy
% balance. Identify the terms that have the greatest and least
% contributions to the energy out of the system. How does this affect
% the efficiency of the engine?
%
% About 85% of losses were miscellaneous losses, 5-6% were frictional
% losses, and the thermal efficiency (work) was only about 9-10%. It
% seems like the miscellaneous losses are very high.
%
% 2b. Write a statement that describes the comparison of the calculated
% efficiency with that of an ideal Otto cycle [quantify the
% discrepancy with a percentage]. What assumptions have you made in
% the air standard analysis that might cause such discrepancies?
%
% I did not calculate the theoretical values. The assumptions that
% would have resulted in discrepencies include perfect isentropic
% compression and expansion and no losses, along with perfect mixing.
%
% 2c. What would you recommend in order to improve the efficiency of the
% engine? For example, how would turbo-charging, inner cooling,
% split-fire spark plugs, and other hardware alterations change the
% efficiency of the engine?
%
% Inner Cooling would increase the efficiency by inreasing the air
% intake charge density. Turbocharging would increase the efficiency
% by recycling heat energy.
%
% 2d. Combustion of fossil fuels releases carbon dioxide, which is a green-
% house gas, into the environment. A diagram of carbon lifecycle is
% shown, illustrating how auto emissions tend to alter the natural
% balance. Suggest 2 solutions for reducing greenhouse gas emissions
% from combustion engines.
%
% 1. The obvious solution is increasing efficiency by minimizing
% losses, so that less fuel is required for equal output
% 2. Another solution is using more effective fuels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%