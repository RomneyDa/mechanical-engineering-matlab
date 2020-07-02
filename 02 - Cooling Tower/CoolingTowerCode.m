% Dallin Romney 
% Cooling Tower Lab

clear, clc, close all;

%%%%%%%%%%%% Data %%%%%%%%%%%%%

    % Constants
    Tamb = 22.0;  % Ambient temperature in lab, Celsius
    Pamb = 656.4*0.133322; % Barometric pressure in lab, kPa
    Qin  = 1.5;   % Input power to water heaters, kW
    Pumpin = 0.1; % Input power to pump, kW
    
    % Small/big t data points at m = 50 g/s
    T50 = [24.85, 15.50, 24.80, 23.10, 27.40, 21.00];
    t50 = [24.30, 23.30, 25.40, 22.00, 21.70, 23.30, 20.45, 21.45, 21.15];

    % 30 g/s
    T30 = [25.20, 15.40, 24.40, 22.40, 33.00, 19.60];
    t30 = [23.60, 22.10, 25.90, 21.40, 21.20, 22.30, 19.40, 22.00, 19.20]; 
    
    % 20 g/s
    T20 = [25.40, 15.30, 24.20, 22.00, 37.40, 19.00];
    t20 = [23.60, 21.40, 25.75, 20.70, 20.50, 21.60, 18.70, 22.30, 18.15];
    
    % Reorganized for for loops:
    T1 = [T50(1), T30(1), T20(1)];
    T2 = [T50(2), T30(2), T20(2)];
    T3 = [T50(3), T30(3), T20(3)];
    T4 = [T50(4), T30(4), T20(4)];
    T5 = [T50(5), T30(5), T20(5)];
    T6 = [T50(6), T30(6), T20(6)];
    
    % Pressures (all trials had same results) [mm H20]
    dPA = 20.5; % Air pressure drop at A
    dPB = 10.0; % Air pressure drop at B
    
    PA = Pamb + 0.101972*dPA; % Pressures at A and B (kPa)
    PB = Pamb + 0.101972*dPB;
        
%%%%%%%%%%%%% Data Reduction %%%%%%%%%%%%%

fprintf('1. Determining dry air mass flow rate:\n\t');


% Psychrometric chart properties [A, B] or [inlet, outlet]
% Specific volume (v, m^3/kg), relative humidity (p), and spec humidity (w)
    v50 = [1.025, 1.045]; p50 = [0.40, 0.88]; w50 = [0.0091, 0.0208];
    v30 = [1.025, 1.041]; p30 = [0.39, 0.84]; w30 = [0.0091, 0.0196];
    v20 = [1.026, 1.040]; p20 = [0.37, 0.83]; w20 = [0.0090, 0.0191];

    disp('1a. See code for data from tables')
    
% Saturation Pressure

    PgIn  = [0, 0, 0]; % Inlets [50 g/s, 30 g/s, 20 g/s]
    PgOut = [0, 0, 0]; % Outlets
    
    for i = 1:3
        Tin  = T1(i) + 273.15; % Convert celsius to Kelvin
        Tout = T3(i) + 273.15;
        
        zIn  = 1 - Tin /647.096;
        zOut = 1 - Tout/647.096;
        
        PgIn(i)  = 22064*exp(647.096/Tin* (-7.8595*zIn  + 1.8441*zIn ^1.5 - 11.7866*zIn ^3 + 22.6807*zIn ^3.5 - 15.9619*zIn^4  + 1.8012*zIn ^7.5));
        PgOut(i) = 22064*exp(647.096/Tout*(-7.8595*zOut + 1.8441*zOut^1.5 - 11.7866*zOut^3 + 22.6807*zOut^3.5 - 15.9619*zOut^4 + 1.8012*zOut^7.5));
    end
    
    fprintf('\t1b. Saturation Pressure [kPa] (50, 30, 20 g/s):\n');
    fprintf('\t\tin: '); disp(PgIn);
    fprintf('\t\tOut: '); disp(PgOut);

% Absolute humidity

    wA = [0, 0, 0];
    wB = [0, 0, 0];
    
    pIn  = [p50(1), p30(1), p20(1)];
    pOut = [p50(2), p30(2), p20(2)];
    
    for i = 1:3
        wA(i) = 0.622*pIn(i) *PgIn (i)/(PA - pIn (i)*PgIn (i));
        wB(i) = 0.622*pOut(i)*PgOut(i)/(PB - pOut(i)*PgOut(i));
    end
    
    fprintf('\t1c. Absolute humidity [kg H20/kg Air] (50, 30, 20 g/s):\n');
    fprintf('\t\tin: '); disp(wA);
    fprintf('\t\tOut: '); disp(wB);
    
% Amagat's Law --> air / water vapor specific volume
    
    Ra = 0.287 ; % Air gas constant (kPa m^3 / kg K)
    Rv = 0.4615; % Water vapor gas constant, same units
    
    vaA = Ra*T1/PA;
    vaB = Ra*T3/PB;
    vvA = Rv*T1/PA;
    vvB = Rv*T3/PB;
    
    vA = vaA + wA.*vvA;
    vB = vaB + wB.*vvB;
    
    fprintf('\t1d. Specific volume [m^3/kg] (50, 30, 20 g/s):\n');
    fprintf('\t\tin: '); disp(vA);
    fprintf('\t\tOut: '); disp(vB);
    fprintf('\t\tThese are all slightly lower than the chart value\n');
    
% Mass flow rate of dry air (from orifice calibration eqn)

    mA = 0.0137*sqrt(dPB./((1 + wB).*vB)')';
    fprintf('\t1f. Mass flow rate of dry air [kg/s] (50, 30, 20 g/s):\n');
    fprintf('\t\t '); disp(mA);    
    
% Specific enthalpy of liquid water [50, 30, 20]

    hIn =  [115.5, 138.3, 156.5]; % kJ/kg
    hOut = [88.1 ,  80.1,  79.8];
    
    fprintf('2. See attached code for table enthalpy values.\n');

% Absolute humidity for state points [ F, G, H ]
    
    wFGH50 = [0.0177, 0.0200, 0.0220];
    wFGH30 = [0.0160, 0.0190, 0.0210];
    wFGH20 = [0.0151, 0.0183, 0.0205];
    
    fprintf('3. See attached code for table humidity values.\n');

% Mass flow rate of water vapor 

    mvA = wA.*mA;
    mvB = wB.*mA;
    
    fprintf('4. Mass flow rates of water vapor [kg/s] (50, 30, 20 g/s):\n');
    fprintf('\tat A: '); disp(mvA);
    fprintf('\tat B: '); disp(mvB);    
    
% Enthalphy of air/water vapor mixture at states A and B

    cP = 1.005; % specific heat of air, 1.005 kJ/kg K

    hgA = 2500.9 + 1.82*T1;
    hgB = 2500.9 + 1.82*T3;
    
    hA = cP*T1 + wA.*hgA;
    hB = cP*T3 + wB.*hgB;
    
    % From chart: (50, 30, 20)
    hAChart = [];
    hBChart = [];
    
    fprintf('5. Enthalpy of Air/vapor mixture [kJ/kg] (50, 30, 20 g/s):\n');
    fprintf('\t\tat A: '); disp(hB);
    fprintf('\t\tat B: '); disp(hB);
    fprintf('\t\tThese are all slightly lower than the chart values\n');
    
    
% All water exits as liquid water --> 

    mIn = [0.05, 0.03, 0.02];
    m6 = mIn + mvA - mvB;
    
% Makeup flow rate

    mMakeup = mIn - m6;

% Total heat gain

    Qair = mA.*(hB - hA);
    Qamb = Qin - mA.*(hB - hA);
   
% Range and approach

%%%%%%%%%%%%%% PRESENTATION %%%%%%%%%%%%%

% Humidity vs height plot
    height = [0, 24.8, 48.3, 71.8]; % cm, (A, F, G, H)

    figure; hold on
    plot(height, [wA(1), wFGH50], 'o-');
    plot(height, [wA(2), wFGH30], '*-');
    plot(height, [wA(3), wFGH20], 'd-');
    
    title('Absolute Humidity vs. Height');
    legend('50 g/s', '30 g/s', '20 g/s');
    ylabel('Absolute Humidity (kg water/kg air)');
    xlabel('Height (cm)');
    
% Dry bulb temp vs Cooling tower height plot

    figure; hold on
    plot(height, [T50(1), t50(8), t50(5), t50(2)], 'o-');
    plot(height, [T30(1), t30(8), t30(5), t30(2)], '*-');
    plot(height, [T20(1), t20(8), t20(5), t20(2)], 'd-');
    
    title('Dry bulb Temp vs. Height');
    legend('50 g/s', '30 g/s', '20 g/s');
    xlabel('Height (cm)');
    ylabel('Dry Bulb Temperature (Celsius)');
    
% m6/m5 vs T5

    figure; plot(T5, m6/mIn, 'r*-')
    title('m6/m5 vs T5');
    xlabel('m6/m5');
    ylabel('T5 (Celsius)');
    
% Qair and Qamb vs. T5 

    figure; hold on;
    plot(T5, Qair, '*-');
    plot(T5, Qamb, 'd-');
    title('Qair and Qamb vs. T5');
    legend('Qair', 'Qamb');
    xlabel('Q (kJ)');
    ylabel('T5 (Celsius)');
    