%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% Fully Developed Pipe Flow                                 %                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%% Experiment properties %%%%%%%%%%%%%%%%%%%%%%%%%%

    Tw = 15;       % Water temperature, degrees Celsius
    mu = 0.0011373;% Dynamic viscosity of water, Pa*s (Engineering Toolbox)
    rho = 999.06;  % Density of water, kg/m^3 (from Engineering Toolbox)
    g = 9.803;     % Acceleration due to gravity, m/s^2
        
    data = xlsread('PipeFlow_ExampleDataSet.xlsx'); % Read given data
    
    % Extract data. All data will be in the form of:
    %         row 1 = 1" PVC, row 2 = 3/4" PVC, row 3 = 3/4" Galvanized
    
    D = [1.033; 0.81; 0.824]*0.0254;                        % Pipe diam, m
    Q = [data(8,2:4); data(8,5:7); data(8,8:10)]*(0.0254^3);% Flow rate m^3
    
    Ldat = data(1:3, 3:7)*0.0254;                           % L1-L5, m
    
    % Populate a 3D array of lengths for simpler calculations later
    L = zeros(3, 3, 5);    % Initate 3D array
    for k = 1:3
        L(:, k, :) = Ldat; % Add layer to 3D array
    end
    
    Ps = zeros(3, 3, 5);   % Same for pressure; 3 pipes, 3 Re's, 5 taps
    for i = 1:3
        for j = 1:3
            Ps(i, j, :) = data(10:14, 3*(i-1) + j + 1);
        end
    end
    Ps  = Ps*0.0254*rho*g; % Convert Ps from in H2O to Pa
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of pipe areas for each flow rate
D = [D'; D'; D']';
A = pi*D.^2/4;    

% Calculate average flow velocity
Vavg = Q./A;

% Reynolds number
ReExp = Vavg.*D*(rho/mu);

% Estimate friction factor between points 3 and 5 for each flow
fExp = 2*(Ps(:,:,3) - Ps(:,:,5)).*D./(rho*Vavg.^2.*(L(:,:,5) - L(:,:,3)));

% Estimate minor loss coefficient and equivalent length of union (2 to 3)
% Assuming negligble length (non-negligble length method commented out)
P23  = Ps(:, :, 2) - Ps(:, :, 3); % Pressure drop from tap 2 to tap 3
L32  =  L(:, :, 3) -  L(:, :, 2); % Length from tap 2 to tap 3

Kt   = 2*P23./(rho*Vavg.^2); % - fExp.*L32./D; % Minor loss coeff
LeD  = 2*P23./(rho*fExp.*Vavg.^2); % - L32./D; % Equivalent Length

meanKt = mean(Kt, 2); % Average union minor loss coefficient for each pipe

% Calculate e/D for each flow for each pipe (solved from Colebrook Eq.)
eDFun = @(f, Re) 3.7*(10.^(-1./(2*sqrt(f))) - 2.51./(Re.*sqrt(f)));

eDExp = eDFun(fExp, ReExp); % e/D for every trial
meaneD = mean(eDExp, 2);    % Average per pipe

%%%%%%%%% Plot of Static Gauge Pressure vs Length/Diameter Ratio %%%%%%%%%%

LD = L./D; % Non-dimensionalized Length for every single trial

% Rotate arrays for plotting purposes - put last dimension first
LDx = permute(LD, [3,2,1]);
Psy = permute(Ps, [3,2,1])/(0.0254*rho*g);

% Array of pipe names and line styles for populating legend
pipeNames = {'1" PVC', '3/4" PVC', '3/4" Galvanized Steel'};
linestyle = {'-ko', '-k*', '-kx', '-go', '-g*', '-gx', '-bo', '-b*','-bx'};
Legend = cell(1, 9); % Initiate empty legend
figure; hold on

for pipe = 1:3 % For each pipe
  for RE = 1:3 % For each Reynolds number
    % Add the current line to the plot and an entry to the legend.
    plot(LDx(:, RE, pipe), Psy(:, RE, pipe), linestyle{3*(pipe-1) + RE});
    Legend{3*(pipe-1) + RE} = [pipeNames{pipe}, ', ', ...
          sprintf('Re = %.2s', ReExp(pipe, RE))];
  end
end

% Annotate
ylim([0, 55]);
legend(Legend);
xlabel('L/D (Dimensionless)');
ylabel('Static Gauge Pressure (in H2O)');

%%%%%%%%%%%%%% Plot of Friction Factor vs. Reynolds Number %%%%%%%%%%%%%%%%

Re = [4e3:1e3:9e3,1e4:1e4:9e4,1e5:1e5:1e6]'; % Reynolds numbers for plot
f = zeros(length(Re), 1);                    % Empty friction factor array

% Colebrook equation, all moved to one side for fzero approximation
colebrook = @(f, Re, eD) 1/sqrt(f) + 2*log10(eD/3.7 + 2.51/(Re*sqrt(f)));

Legend = cell(1, 12);           % Initiate empty legend
pipeStyle = {'ko', 'gs', 'bd'}; % Line styles for plotting in the loop

figure; 
for j = 1:length(meaneD)  % For each e/D surface roughness value
    for k = 1:length(Re)  % For each reynolds number
        % Calculate the friction factor by zeroing the above equation
        f(k,1) = fzero(@(f) colebrook(f,Re(k),meaneD(j)),0.026);
    end
    loglog(Re,f,pipeStyle{j}(1));    % Add to diagram
    hold on
    Legend{j} = pipeNames{j};        % Add entry to legend
end
grid on

meanReExp = mean(ReExp, 2)';     % Desired Reynolds number values

% Annotate and adjust y limits
xlabel('Reynolds Number Re');
ylabel('Friction Factor f');
ylim([0.015, 0.05]);

% Add data points onto plot
for pipe = 1:3 % For each pipe
  for RE = 1:3 % For each Reynolds number
    % Add the current point to the plot
    loglog(ReExp(pipe, RE), fExp(pipe, RE), pipeStyle{pipe});
    
    % Add a legend entry
    Legend{3*(pipe-1) + RE + 3} = [pipeNames{pipe}, ', ', ...
         sprintf('Re = %.2s', ReExp(pipe, RE))];
  end
end

% Annotate
legend(Legend);

%%%%%%%%%%%%%%%%%%%% Other calculations for discussion %%%%%%%%%%%%%%%%%%%%

% Determine theoretical union K values from empirical correlation 
KCouplings = 0.083*(mean(D/0.0254, 2)).^(-0.69);

eValues    = 100*mean(eDExp, 2).*mean(D, 2); % Avg Roughness lengths, cm
eVal = 100*eDExp.*D;                         % (For every trial)

LMajor = L(:, :, 5) - L32;                    % Major L (excludes union)
fMajor = fExp.*rho.*(Vavg.^2).*LMajor./(2*D); % Major friction losses
fMinor = fExp.*rho.*(Vavg.^2).*LeD;           % Minor friction losses

% Power required at highest flow rate
Power = Q.*Ps(:, :, end);