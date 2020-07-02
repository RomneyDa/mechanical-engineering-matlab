%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% Transient Conduction                                       %                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

% Material Properties:
cpB = 377;   % Specific heat of brass, J/kg*K
cpA = 875;   % Specific heat of aluminum, J/kg*K
dB  = 8498;  % Density of 360 brass, kg/m^3
dA  = 2780;  % Density of aluminum, kg/m^3
kB  = 116.0; % Thermal conductivity of brass, W/m*K
kA  = 121.4; % Thermal conductivity of aluminum, W/m*K

alphaA = kA/(dA*cpA); % Thermal diffusivity of aluminum (W/m^2 K)
alphaB = kB/(dB*cpB); % Thermal diffusivity of brass (W/m^2 K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Plot 1: Measured Temperatures %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the Transient data for each material
    Brass = xlsread('data/Brass.xls');       
    Aluminum = xlsread('data/Aluminum.xls'); 

% Separate columns of time and temperatures for each material
    % Brass
    BBathTemp  = Brass(:, 1); % Bath temperature (deg C)
    BShapeTemp = Brass(:, 2); % Sphere internal temperature (deg C)
    BTime      = Brass(:, 3); % Times (seconds)     

    % Aluminum
    ABathTemp  = Aluminum(:, 1); % Bath temperature (deg C)
    AShapeTemp = Aluminum(:, 2); % Sphere internal temperature (deg C)
    ATime      = Aluminum(:, 3); % Times (seconds)

% Plot of temperature vs. time for each material
    figure; hold on
    plot(ATime, AShapeTemp, 'o'); % Plot aluminum data
    plot(BTime, BShapeTemp, 'o'); % Plot brass data

    % Add straight lines at the average bath temperatures for each material
    TinfA = mean(ABathTemp);
    TinfB = mean(BBathTemp);

    line([min(ATime), max(ATime)],[TinfA, TinfA],'LineStyle','--');
    line([min(BTime), max(BTime)],[TinfB, TinfB],'Color','r',...
         'LineStyle','--');

    % Annotate
    xlabel('Time (sec)');
    ylabel('Temperature (degrees Celsius)');
    legend('Aluminum Sphere', 'Brass Sphere', 'Aluminum Bath', 'Brass Bath',...
           'LOCATION', 'SOUTHEAST');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plots 2 & 3: Non-dimensional temperatures vs. Fourier number %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r0 = 51/2/1000; % Sphere radius, converted from mm diameter to m

% Non-dimensional temperatures = (T0 - Tinf)/(Ti - Tinf)
    NDTBrass = (BShapeTemp - BBathTemp)./(-BBathTemp + BShapeTemp(1));
    NDTAluminum = (AShapeTemp - ABathTemp)./(-ABathTemp + AShapeTemp(1));

% Fourier Numbers for Brass and Aluminum
    FoB = alphaB*BTime/r0^2;
    FoA = alphaA*ATime/r0^2;

% Calculate a 1st order polynomial fit for linearized sphere conduction eq:
% ln(NDT) = -Z1^2*Fo + ln(C1) (y = m*x + b) from NDT = C1*exp(-Z1^2*Fo)
% y = ln(NDT), m = -Z1^2, x = Fo, b = ln(C1)

    % For each material:
    logNDTA = log(NDTAluminum); % Left side of equation: y = ln(NDT)
    logNDTB = log(NDTBrass); 
    
    indexA = ~isinf(logNDTA);   % Remove infinite values
    indexB = ~isinf(logNDTB);
   
    coeffA = polyfit(FoA(indexA), logNDTA(indexA), 1); % 1st order
    coeffB = polyfit(FoB(indexB), logNDTB(indexB), 1); % polynomial fit

    Z1A = sqrt(-coeffA(1)); % m = -Z1^2 --> Z1 = sqrt(-m); 1st coefficient
    Z1B = sqrt(-coeffB(1));

    C1A = exp(coeffA(2));   % b = ln(C1) --> C1 = e^b 1st coeff in fit
    C1B = exp(coeffB(2));
    
    FitAluminum = C1A*exp(-Z1A^2*FoA); % Reinsert into sphere equation
    FitBrass = C1B*exp(-Z1B^2*FoB);    % aka nonlinearize

% Plot and annotate nondimensional time vs. fourier number and accompanying
% fit for brass and aluminum (separately), semilog on the y axis, linear on
% the x axis

    % Aluminum
    figure
    semilogy(FoA, NDTAluminum, 'bo', FoA, FitAluminum,  'r-');
    xlabel('Fourier Number Fo');
    ylabel('Nondimensional Temperature');
    legend('Data', 'First order fit');

    % Brass
    figure
    semilogy(FoB, NDTBrass, 'ro', FoB, FitBrass, 'b-');
    xlabel('Fourier Number Fo');
    ylabel('Nondimensional Temperature');
    legend('Data', 'First order fit');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Short Answer Question Stuff %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goodness of fit: Fo number chosen where plot appears to start following
% exponential curve
    [RA, ~] = corrcoef(NDTAluminum(FoA > 0.5), FitAluminum(FoA > 0.5));
    [RB, ~] = corrcoef(NDTBrass   (FoB > 0.4), FitBrass   (FoB > 0.4));

    fprintf('For aluminum, r^2 = %.4f after Fo = 0.5.\n', RA(2)^2);
    fprintf('For brass, r^2 = %.4f after Fo = 0.4.\n'   , RB(2)^2);

% Biot number from Z1 (derivation on handout)
    BiA = 1 - Z1A*cot(Z1A);
    BiB = 1 - Z1B*cot(Z1B);

% Correlation coefficient
    hA = BiA*kA/r0;
    hB = BiB*kB/r0;

% Determine thermal conductivity by using the convection coefficient for
% the OTHER material: based on knowledge that h should be the same for the 
% same shape and fluid properties.
    kAMeasured = hB*r0/BiA;
    kBMeasured = hA*r0/BiB;

% Verification of lumped capacitance method: 
    % First, plot measured data for aluminum, as on the very first plot
    figure; hold on
    plot(ATime, AShapeTemp, 'o');
    line([min(ATime), max(ATime)],[mean(ABathTemp), mean(ABathTemp)], ...
        'LineStyle', '--');
    xlabel('Time (sec)');
    ylabel('Temperature (degrees Celsius)');

    % Second, calculate and plot lumped capacitance approximation
    T = TinfA + (min(AShapeTemp) - TinfA)*exp(-BiA*FoA);
    plot(ATime, T);
    legend('Data', 'Bath Temp', 'Lumped Capacitance Approximation',...
    'LOCATION','SOUTHEAST');