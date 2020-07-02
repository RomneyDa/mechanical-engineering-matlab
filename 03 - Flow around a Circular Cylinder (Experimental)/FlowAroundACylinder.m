% Dallin Romney
% Flow Around a Circular Cylinder

clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = 0.75*0.0254;     % Cylinder diameter, m
w = 12*0.0254;       % Wind tunnel width, m

Pbar = 29.9; % Barometric pressure (inHg converted to atm)
T = 21.3 + 273.15;   % Air temperature (K)
R = 0.00008206;      % Universal gas constant (m^3*atm/(mol*K))
m = 0.02897;         % Molar mass of air, kg/mol
nu = 0.00001523;     % Kinematic viscosity of air (m^2/s)

rho = (Pbar*0.03342)*m/(R*T);  % Air density (kg/m^3)

Pdynam  = mean(textread('data/Approach_Dyn.txt'));    % Mean measured dynamic pressure
Pstatic = mean(textread('data/Approach_Static.txt')); % Mean measured static pressure

C = 4.7538; % Bernoulli’s equation const. (m*s^{-1}*in^{1/2}*(mm*K)^{-1/2})

heights = [17.79, 16.76, 15.76, 14.78, 13.75, 13.50, 13.25, 13.00, ...
           12.74, 12.50, 12.25, 12.00, 11.75, 11.50, 11.25, 11.00, ...
           10.00, 9.00, 8.00, 7.00]; % Vertical pos of measurement, cm

measureDist = 25;                    % Horizontal pos from cylinder, cm
vectorLength = 602;                  % Number of measurements/trial
numFiles = length(heights);          % Number of 30 second trials

pressure = zeros(vectorLength, numFiles); % Initiate array for data

for k = 1:numFiles                                % For each cile
    filename = sprintf('data/Wake_Dynamic_%d.txt', k); % Current file name 
    pressure(:, k) = textread(filename, '%f');    % Scan file into column
end

meanPressure = mean(pressure); % Average velocities for each file

[~, midPos] = min(meanPressure); % Find midpoint by finding minimum of vel
midHeight = heights(midPos);     % Height that matches up with cylinder

%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT of u vs. y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = C*sqrt(meanPressure*T/Pbar); % Calculate velocity
y = (heights - midHeight)/100;   % Center the height values, convert to m

plot(u, y); % Plot velocity profile (m/s vs m)

tit = sprintf('Vertical Position in Wake vs. Velocity, %d cm downstream', measureDist);
title(tit); 
ylabel('Vertical Distance from Center of Wake (m)'); 
xlabel('Air Velocity (m/s)'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT of Re vs. Cd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velocityFreeStream = sqrt(Pdynam*133.32*2/rho); % Free stream velocity
Re = velocityFreeStream*D/nu ;                  % Reynolds number around cylinder

[Re_book, Cd_book] = textread('data/CD_RE_Textbook.txt','%f %f');

figure
semilogx(Re_book, Cd_book);

title('Reynolds Number for Cylinder vs. Drag Coefficient'); 
xlabel('Reynolds Number'); 
ylabel('Drag Coefficient'); 
legend('Textbook');
