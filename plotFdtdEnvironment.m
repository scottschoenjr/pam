%**************************************************************************
%
% FDTD Environment Plotter
%
%   Function plots data created with the FDTD data generator. Loads in data
%   plots the sound speed, density, and attenuation fields defined therein.
%
%              Scott Schoen Jr | Georgia Tech | 20170105
%
%************************************************************************** 

function [] = plotFdtdEnvironment( dataFile ) 

% Simulation parameters
recPosition = 100E-3; % Distance from x = 0 [m]

% Define sources
% sourcePositions = [ ... % (x, y) positions [m]
%     40, 40; ...
%     40, 41.8; ...
%     40, 38.2 ...
%     ]./1E3;
sourcePositions = [ ... % (x, y) positions [m]
    40, 40 ...
    ]./1E3;

% Load in the data file of interest
data = load(dataFile);
domain = data.domain;

% Get relevant parameters from file
dx = domain.res*1E-3; % Distance between nodes [m]
rho = flip( flip(domain.extrhotot, 2)', 2);   % Density
c = flip( flip(domain.extctot, 2)', 2);       % Sound speed
alpha = flip( flip(domain.extatttot, 2)', 2); % Attenuation

% Get dimensions and position vectors
xDim = size(rho, 2); % Number of nodes in x direction
yDim = size(rho, 1); % Number of nodes in y direction

xPositionVector = (0:dx:(xDim-1)*dx); % [m]
yPositionVector = (0:dx:(yDim-1)*dx); % [m]

[x, y] = meshgrid( xPositionVector, yPositionVector );


figure()
hold all;

% Plot the sound speed field
pcolor( x.*1E3, y.*1E3, c );
shading flat;

% Plot the sources' and receiver positions
% plot( sourcePositions(:, 1).*1E3, sourcePositions(:, 2).*1E3, 'ro' );
plot( 1E3.*[recPosition, recPosition], ...
    1E3.*[min(yPositionVector), max(yPositionVector)], '--w' );

% Formatting
% xlabel( 'Distance from Receiver [mm]' );
% ylabel( 'Sensor Position [mm]' );
% 
% xlim( 1E3.*[min(xPositionVector), max(xPositionVector)] );
% ylim( 1E3.*[min(yPositionVector), max(yPositionVector)] );
% 
% cBarHandle = colorbar;
% 
% cBarHandle.Label.String = 'Sound Speed [m/s]';
% cBarHandle.Label.FontSize = 16;
% cBarHandle.Label.Interpreter = 'latex';
% cBarHandle.TickLabelInterpreter = 'latex';

xlabel( 'Axial Distance [mm]', 'FontSize', 26 );
ylabel( 'Transverse Distance [mm]', 'FontSize', 26 );

xlim( 1E3.*[min(xPositionVector), max(xPositionVector)] );
ylim( 1E3.*[min(yPositionVector), max(yPositionVector)] );
xlim([0, 120]);
axis equal;

cBarHandle = colorbar;

cBarHandle.Label.String = 'Sound Speed [m/s]';
cBarHandle.Label.FontSize = 26;
cBarHandle.Label.Interpreter = 'latex';
cBarHandle.TickLabelInterpreter = 'latex';


end