%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% Contraflow Heat Exchangers                                %                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%% Experiment properties %%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
    
    Qh   = [3.0, 2.5, 1.5, 3.0]/15850.3; % Hot  volume flow rate, m^3/s
    Qc   = [5.0, 5.0, 5.0, 3.0]/15850.3; % Cold volume flow rate, m^3/s
    
    Thi  = ([115.3,114.9,115.1,116.9] - 32)*5/9 + 273.15; % Hot temp in, K
    Tho  = ([100.4, 97.7, 91.3,105.1] - 32)*5/9 + 273.15; % Hot  out, K
    Tci  = ([ 54.6, 53.8, 53.8, 53.9] - 32)*5/9 + 273.15; % Cold  in, K
    Tco  = ([ 64.1, 63.0, 62.0, 66.4] - 32)*5/9 + 273.15; % Cold out, K
    
    Tsh  = ([103.3,  95.7,  94.4,  97.2] - 32)*5/9 + 273.15; % Shell temp
    Tamb = ([ 70.2,  70.6,  70.6,  70.6] - 32)*5/9 + 273.15; % Ambient temp
    
    % The specific heat is very constant, varies from about 4180 to 4190  
    % j/kgK over the entire temperature range. Density ranges from about 
    % 989 kg/m^3 to 999 kg/m^3, so an average is taken. The maximum error 
    % is 0.5 %. This is plenty good for these calculations.    
    rho = 994;  % Density of water, kg/m^3 at average temperature
    cp  = 4186; % Isobaric specific heat of water, J/kgK
    
    % Heat exchanger properties
    F  = 1;  % Correction factor, counterflow single pass shell and tube
    n  = 31; % Number of tubes
    
    Ds = 2.120*0.0254; % Shell diameter, m
    L  = 9.000*0.0254; % Tube length, m
    bs = 1.130*0.0254; % Baffle spacing, m
    wt = 0.028*0.0254; % Tube wall thickness, m
    Do = 0.250*0.0254; % Tube outer diameter, m
    
    % Heat transfer area, m^2, calculated at avg diameter
    A0 = pi*(Do - wt/2)*L*n;
    
    % Radiation and convection calculations
    emiss = 0.95; % Emissivity (painted gray)
    As = pi*Ds*L; % Shell area, m^2
    
    Tf  = (mean(Tamb) + mean(Tsh)) / 2; % Film temperature
    kf  = 0.0257;    % Thermal conductivity of air at 294.5 K, W/mK
    Pr  = 0.71393;   % Prandtl number at film temperature
    v   = 0.0000159; % Kinematic visocity of air at film temperature
    V   = 0.5;       % Airflow speed, m/s. Assumption made for convection
    ReD = V*Ds/v;    % Reynolds number of flow
    
    % Nusselt number for horizontal flow over a circular cylinder
    NuD = 0.3 + (0.62*ReD^0.5*Pr^(1/3)/(1 + (0.4*Pr)^(2/3))^(1/4)) * ...
                (1 + (ReD / 282000)^(5/8))^(4/5); 
    
% Data Reduction
nt = length(Qh); % number of trials

% Initiate vectors for for loop
q    = zeros(1, nt); % Heat transfer rate on cold side
qh   = zeros(1, nt); % Heat transfer rate on hot side
dq   = zeros(1, nt); % Percent difference in heat transfer rates
Cc   = zeros(1, nt); % Cold heat capacity rate
Ch   = zeros(1, nt); % Hot  heat capacity rate
Cmin = zeros(1, nt); % Smaller heat capacity rate
Cmax = zeros(1, nt); % Larger  heat capacity rate
Cr   = zeros(1, nt); % Heat capacity ratio
Tlm  = zeros(1, nt); % Log mean temp difference
U0   = zeros(1, nt); % Overall heat transfer coefficient
NTU  = zeros(1, nt); % Number of transfer units
effM = zeros(1, nt); % Measured effectiveness
effT = zeros(1, nt); % Theoretical effectiveness
deff = zeros(1, nt); % Percent difference in measured/theo effectivenesses
ThoT = zeros(1, nt); % Theoretical hot temp out
dTho = zeros(1, nt); % Percent difference in measured/theo hot temp out
qrad = zeros(1, nt); % Radiative heat transfer
qcnv = zeros(1, nt); % Convective heat transfer

for trial = 1:nt     % For each trial,
        
    % Calculate each of the values listed above
    Cc(trial) = Qc(trial)*rho*cp;
    Ch(trial) = Qh(trial)*rho*cp;
     q(trial) = Cc(trial)*(Tco(trial) - Tci(trial));
    qh(trial) = Ch(trial)*(Thi(trial) - Tho(trial));
    dq(trial) = 100*abs(qh(trial) - q(trial)) / q(trial);
     
  Cmin(trial) = min([Cc(trial), Ch(trial)]);
  Cmax(trial) = max([Cc(trial), Ch(trial)]);
    Cr(trial) = Cmin(trial) / Cmax(trial);
    
          dt1 = Thi(trial) - Tco(trial);
          dt2 = Tho(trial) - Tci(trial);
   Tlm(trial) = (dt2 - dt1) / log(dt2/dt1);
    U0(trial) = q(trial) / (A0*F*Tlm(trial));
   NTU(trial) = U0(trial)*A0 / Cmin(trial);
   
   if Cr(trial) == 1
       effT(trial) = NTU(trial) / (1 + NTU(trial));
   else
       effT(trial) = (1 -           exp(-NTU(trial)*(1 - Cr(trial)))) / ...
                     (1 - Cr(trial)*exp(-NTU(trial)*(1 - Cr(trial))));
   end
   
   effM(trial) = q(trial) / (Cmin(trial)*(Thi(trial) - Tci(trial)));
   deff(trial) = 100*abs(effM(trial) - effT(trial)) / effT(trial);
   
   ThoT(trial) = Thi(trial)*(1 - effM(trial)*Cmin(trial) / Ch(trial)) + ...
                 Tci(trial)*(    effM(trial)*Cmin(trial) / Ch(trial));
   dTho(trial) = 100*abs(ThoT(trial) - Tho(trial)) / Tho(trial);
   
   % Convection
   h = NuD*kf / Ds;
   qcnv(trial) = h*As*(Tsh(trial) - Tamb(trial));
   
   % Radiation
   qrad(trial) = emiss*5.67E-8*(Tsh(trial)^4 - Tamb(trial)^4)*As;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%% Plot of Effectiveness vs NTU %%%%%%%%%%%%%%%%%%%%
openfig('EffectivenessNTU.fig');
figure(1)
hold on

% Overlay experimental values
plot(NTU(1), effM(1), 'bo');
plot(NTU(2), effM(2), 'bo');
plot(NTU(3), effM(3), 'ms');
plot(NTU(4), effM(4), 'ms');

legend('Theoretical', 'Theoretical', 'Theoretical', 'Theoretical', ...
       'Theoretical', 'Trial 1', 'Trial 2', 'Trial 3',  'Trial 4', ...
       'location', 'southeast'); 

