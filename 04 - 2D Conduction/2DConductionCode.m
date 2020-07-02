%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dallin Romney                                                           %
% 2D Conduction                                              %
% October 19, 2018                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

plateSize = 15.24; % length of side of plate, cm (it's a square)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Plot 1: Measured Data Filled Contours %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = xlsread('Section15.csv');

[rows, cols] = size(data);        % Get size of data
x = linspace(0, plateSize, cols); % Create vector of x position values
y = linspace(0, plateSize, rows); % Create vector of y position values
[xMesh, yMesh] = meshgrid(x, y);  % 2D array of x values & same for y vals

% Generate contour plot. The contourf function flips the data vertically
% about the origin so I flipped it back for visual purposes. The "100"
% denotes 100 different colors used, and contour lines are removed.

contourf(xMesh, yMesh, flipud(data), 100, 'LineStyle', 'none'); 
cb = colorbar; cb.Label.String = 'Temperature (deg C)';         % Colorbar
xlabel('X (cm)'); ylabel('Y (cm)');                             % Annotate
%title('Measured');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Plot 2: Computed Values - Adiabatic Bottom Filled Contours %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate plate with boundary conditions. Bottom boundary is adiabatic, so
% it will stay at zeros. Right, left and top sides (and all the corners)
% are replaced with measured values.

T = zeros(rows, cols);    % Matrix of Temperatures
T(:,   1) = data(:,   1); % Replace left side with measured boundary 
T(:, end) = data(:, end); % Replace right side with measured boundary 
T(1,   :) = data(1,   :); % Replace top row with measured boundary

tol = 0.000002; % Tolerance of Gauss-Seidel Method
rel = 2*tol;    % Initiate the relative difference to greater than tol

% While the maximum relative difference is greater than the tolerance
while (max(max(rel)) > tol)  
    
    Told = T;
    
    % Determine the new temperature at the interior nodes
    for j = 2:cols - 1
       for i = 2:rows - 1
          T(i, j) = (T(i, j+1) + T(i, j-1) + T(i-1, j) + T(i+1, j)) / 4;
       end
    end
    
    % Determine the new temperature at the adiabatic boundary (bottom)
    for j = 2:cols - 1
        T(rows, j) = (T(rows, j+1) + T(rows, j-1) + 2*T(rows-1, j)) / 4;
    end
    
    % Calculate the relative difference between iterations (rel = array)
    rel = (T - Told) ./ Told;
    
end

% Plot the calculated values in the same manner as the measured ones
figure
contourf(xMesh, yMesh, flipud(T), 100, 'LineStyle', 'none'); % Contour plot
cb = colorbar; cb.Label.String = 'Temperature (deg C)';      % Add colorbar
xlabel('X (cm)'); ylabel('Y (cm)');                          % Annotate
%title('Computed: Adiabatic Boundary Condition');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Plot 3: Measured Data Isolines %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numLines = 17;                            % Number of isolines

maxT = max(max(data));                    % Maximum measured temperature
minT = min(min(data));                    % Minimum measured temperature    
isoVals = linspace(minT, maxT, numLines); % Values at which to draw lines

figure
[C, hnd] = contour(xMesh, yMesh, data, isoVals); % Plot isotherms

set(hnd,'color','k','linestyle','-');            % Set isoline linestyle
clabel(C, hnd);                                  % Turn on labels
xlabel('X (cm)'); ylabel('Y (cm)');
%title('Measured Isotherms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Plot 4: Computed Data Isolines %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = max(max(T));                       % Maximum computed temperature
minT = min(min(T));                       % Minimum computed temperature    
isoVals = linspace(minT, maxT, numLines); % Values at which to draw lines

figure
[C, hnd] = contour(xMesh, yMesh, T, isoVals); % Plot isotherms

set(hnd,'color','k','linestyle','-');         % Set isoline linestyle
clabel(C, hnd);                               % Turn on labels
xlabel('X (cm)'); ylabel('Y (cm)');
%title('Computed Isotherms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot 5: Difference between Computed and Measured - Filled Contours %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
contourf(xMesh, yMesh, flipud(data - T), 100, 'LineStyle', 'none');
cb = colorbar; cb.Label.String = 'Temperature Difference (deg C)';
xlabel('X (cm)'); ylabel('Y (cm)');
%title('Difference between Measured and Computed (M-C)');




