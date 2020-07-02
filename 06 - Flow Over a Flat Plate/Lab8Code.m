%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% Convection over a Flat Plate                               %                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

% Experiment properties
    P = 657*133.3; % Ambient pressure, mmHg to Pa
    R = 156.5;     % Resistance, Ohms
    V = 35.89;     % Voltage, V
    f = 9;         % Wind tunnel fan frequency, Hz
    
    Tamb = 21.0 + 273.15; % Ambient Temperature, K

    plateW =  68/1000;      % Width of heated plate, m
    xStart =  77/1000;      % x location of beginning of heated plate, m
    xEnd   = 230/1000;      % x location of end of heated plate, m
    plateL = xEnd - xStart; % Length of heated plate, m

    % x locations of each thermocouple, m
    x = [85 92 102 112 123 134 143 153 162 173 186 196 209 219]/1000;
    
    u = 0.704*f - 1.373; % Free-stream velocity, m/s
    Q = V^2/R;           % Power into plate (two-sided), W
    A = plateW*plateL;   % Plate surface area, m^2
    q = Q/(2*A);         % Heat flux through plate (2 sided), W/m^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Import Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data = textread('Flat Plate.txt', '%f'); % Read in measurements
    T = data(2:2:end)' + 273.15;             % Extract temperatures (K)
    T6 = T(6); T12 = T(12);                  % Extract redundancies
    T(6) = []; T(12) = [];                   % (6 and 12 for conductance)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Properties of Air %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Tf = (T + Tamb)/2; % Film temperature
    
% Emperical / theoretical property correlations: anonymous functions
  % Thermal conductivity of air from empirical correlation
    kFun = @(T) T.^3*1.5207e-11 - T.^2*4.8574e-8 + T*1.0184e-4 - 3.9333e-4;
  % Specific heat of air from empirical correlation
    cFun = @(T) T.^4*1.933e-10 - T.^3*8e-7+T.^2*1.141e-3 - T*4.489e-1+1058;
  % Density of air from ideal gas law
    dFun = @(T, p) p./(287.05*T);
  % Dynamic viscosity of air from the Sutherland equation
    vFun = @(T) (T.^(1.5)*1.458e-6)./(T + 110.4);
    
    % Experimental thermal conductivities and densities/viscosities for air
    kf = kFun(Tf);
    d  = dFun(Tf, P);
    v  = vFun(Tf);
    cp = cFun(Tf);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Experimental values %%%%%%%%%%%%%%%%%%%%%%%%%

% Plot measured surface temperatures
figure; hold on
plot(x*1000, T - 273.15, 'bo');
xlabel('x position, mm'); ylabel('Temperature (deg C)');

% Experimental convection coefficients; plot vs. x
hExp = q./(T - Tamb);
figure; hold on
plot(x*1000, hExp, 'bo');
xlabel('x position, mm'); ylabel('Convection Coefficient (W/m^2 K');

% Experimental Nusselt numbers at each x value; plot vs. x
NuExp = hExp.*x./kf;
figure; hold on
plot(x*1000, NuExp, 'bo');
xlabel('x position, mm'); ylabel('Nusselt Number');

%%%%%%%%%%%%%%%%%%%%%%%%%%% Theoretical values %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prx = cp.*v./kf;  % Prandlt number at each x value
Rex = d.*u.*x./v; % Reynolds number at each x value

% Nusselt number from empirical correlation; add to plot
Nux = (0.453*Rex.^0.5.*Prx.^(1/3))./((1 - (xStart./x).^0.75).^(1/3));
figure(3); plot(x*1000, Nux, 'r-');
legend('Experimental', 'Theoretical', 'Location', 'Northwest');

% hx from empirical nusselt number; add to plot
hx = Nux.*kf./x;
figure(2); plot(x*1000, hx, 'r-');
legend('Experimental', 'Theoretical');

% Theoretical surface temperatures based on measured power and theo h; plot
Tx = q./hx + Tamb;
figure(1); plot(x*1000, Tx - 273.15, 'r-');
legend('Experimental', 'Theoretical', 'Location', 'Southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RADIATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    emiss = 0.7;                             % Plate emissivity
    qRad = emiss*(5.67e-8)*(T.^4 - Tamb.^4); % Radiative heat flux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The widths used for integration are the averages of consecutive
% differences between x positions. The first and last ones use the distance
% from the edge plus half the distance to the closest value.
dx = diff(x);
xW = [x(1)-xStart+dx(1)/2,(dx(1:end-1)+dx(2:end))/2,xEnd-x(end)+dx(end)/2];

% Anonymous functions for averaging and integrating
integrateX = @(var) sum(xW.*var) * plateW; % Int = sum(width*var*dx)
averageX   = @(var) sum(xW.*var) / plateL; % Avg = 1/L*int(var*dx)

% Weighted averages along x direction
    hTheoAvg  = averageX(hx);              % Convection coefficients
    hExpAvg   = averageX(hExp);            %
    kfAvg     = averageX(kf);              % Thermal conductivity
    NuTheoAvg = hTheoAvg*plateL / kfAvg;   % Nusselt numbers ( = f(h, kf))
    NuExpAvg  = hExpAvg *plateL / kfAvg;   %
    TTheoAvg  = averageX(Tx);              % Temperatures
    TExpAvg   = averageX(T);               %
    qRadAvg   = averageX(qRad);            % Radiation flux
    
% Heat loss due to radiation (two sides)
    QRad      = 2*integrateX(qRad);       
    
% Simple arithmetic means (yields similar results to weighted averages)
    NuTheoMean = mean(Nux);                % Nusselt numbers
    NuExpMean  = mean(NuExp);              %
    hTheoMean  = mean(hx);                 % Convection coefficients
    hExpMean   = mean(hExp);               %
    TTheoMean  = mean(Tx);                 % Temperatures
    TExpMean   = mean(T);                  %
    