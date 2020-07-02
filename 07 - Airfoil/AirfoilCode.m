%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% Airfoils                                                   %                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%% Experiment properties %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Air properties
    P = 657.1*133.3;              % Ambient pressure, mmHg to Pa
    Tamb = 21.7 + 273.15;         % Ambient Temperature, K
    
    dFun = @(T, p) p./(287.05*T); % Density empirical correlation
    rho = dFun(Tamb, P);          % Air density
   
    muFun = @(T) (T.^(1.5)*1.458e-6)./(T + 110.4); % Dyn visc (sutherland)
    mu = muFun(Tamb);             % Dynamic viscosity of air
    
    Pdyn  = 1.36*248.8;   % Dynamic pressure, in H20 conv to Pa
    Pstat = 2.46*248.8;   % Static  pressure, in H20 conv to Pa
    
    V = sqrt(2*Pdyn/rho); % Air velocity, m/s
    
    % Data for lift and drag measurements
    drag0 = -0.132; % Initial drag, kg
    angle = [0, 3, 5, 8, 10, 12, 14, 16];
    lift =  [0.136, 0.363, 0.430, 0.516, 0.452 ,0.410, 0.402, 0.412];
    drag = -[0.149, 0.150, 0.150, 0.160, 0.245, 0.287, 0.310, 0.325]-drag0;
    
    % Data for static pressure measurements, pressures are in in H20
    x   = [3, 9.5, 20, 30, 40, 55, 60, 70, 80]/1000; % x position in m
    P_5 =  [0.62, 0.12, -0.11, -0.18, -0.20, -0.21, -0.19, -0.17, -0.13];
    P5  = -[2.19, 1.53,  1.16,  0.84,  0.60,  0.48,  0.36,  0.23,  0.10];
    P15 = -[0.78, 0.80,  0.79,  0.77,  0.78,  0.78,  0.80,  0.82,  0.85];
    
    % Airfoil properties
    s = 12*0.0254; % Span length, in converted to m
    c =  4*0.0254; % Chord length, in converted to m
    Ap = s*c;      % Planform area, m^2
    
    % Reynolds number
    Re = V*c*rho/mu;
    
%%%%%%%%%%%%%%%% Plot of Lift Coefficient vs Attack Angle %%%%%%%%%%%%%%%%%
    
    FL =  lift*9.81; % Lift force, kg converted to N
    FD = -drag*9.81; % Drag force, kg converted to N
 
    CL = 2*FL/(rho*V^2*Ap); % Lift coefficients
    CD = 2*FD/(rho*V^2*Ap); % Drag coefficients
    
    % Plot CL vs angle and annotate
    plot(angle, CL);
    xlim([0, 18]);
    ylim([0, max(CL)]);
    xlabel('Angle of Attack (degrees)');
    ylabel('Lift Coefficient CL');
    
%%%%%%%%%%%%%%%% Plot of Drag Coefficient vs Attack Angle %%%%%%%%%%%%%%%%%
   
    % On a new figure, plot CD vs angle and annotate
    figure
    plot(angle, CD);
    xlim([0, 18]);
    ylim([0, max(CD)]);
    xlabel('Angle of Attack (degrees)');
    ylabel('Drag Coefficient CD');
    
%%%%%%%%%%%%%% Plot of Lift Coefficient vs Drag Coefficient %%%%%%%%%%%%%%%

    figure; hold on
    plot(CD, CL, 'bo'); % Plot experimental CL vs CD from lab
    
    data = xlsread('NACA0012_Data.xlsx'); % Accepted NACA0012 airfoil data

    CL_Re1 = data(1:19, 2); % Lift coefficient at Re = 1*10^6
    CD_Re1 = data(1:19, 3); % Drag coefficient at Re = 1*10^6
    CL_Re2 = data(1:19, 5); % Lift coefficient at Re = 2*10^6
    CD_Re2 = data(1:19, 6); % etc.
    CL_Re5 = data(1:19, 8);
    CD_Re5 = data(1:19, 9);
    
    % Add accepted values at various Reynolds numbers onto plot, annotate
    plot(CD_Re1, CL_Re1, 'r--');
    plot(CD_Re2, CL_Re2, 'k--');
    plot(CD_Re5, CL_Re5, 'g--');
    
    ylabel('Lift Coefficient CL');
    xlabel('Drag Coefficient CD');
    
    message = sprintf('Experimental (Re = %d)', Re);
    legend(message, 'Re = 1*10^6', 'Re = 1*10^6', 'Re = 5*10^6', ...
           'location', 'northeast');
       
%%%%%%%% Plot of Pressure Coefficient vs x/c for angle below stall %%%%%%%%

CP_5 = -2*(P_5*248.8 - Pstat)/(rho*V^2); % Pressure coefficients at -5 deg
CP5 =  -2*( P5*248.8 - Pstat)/(rho*V^2); % Pressure coefficients at  5 deg


% Plot pressure coefficient vs x/c, annotate
figure;  plot(x/c, CP_5);
hold on; plot(x/c, CP5);

message = sprintf('Cp vs x/c for Angle of Attack = 5 deg, Re = %d', Re);
title(message);
xlabel('Ratio of position to chord length x/c');
ylabel('Pressure Coefficient Cp');
legend('Bottom of airfoil', 'Top of airfoil');

%%%%%%%% Plot of Pressure Coefficient vs x/c for angle above stall %%%%%%%%

CP15 = -2*(P15*248.8 - Pstat)/(rho*V^2); % Pressure coefficients at 15 deg

% Plot pressure coefficient vs x/c, annotate
figure;  plot(x/c, CP15);
hold on; plot(x/c, CP_5);

message = sprintf('Cp vs x/c for Angle of Attack = 15 deg, Re = %d', Re);
title(message);
xlabel('Ratio of position to chord length x/c');
ylabel('Pressure Coefficient Cp');
legend('Bottom of airfoil', 'Top of airfoil');

    
