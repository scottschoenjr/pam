%**************************************************************************
%
% Angular Spectrum Field Reconstruction (2D)
%
%   This script uses a passive acoustic mapping (PAM) based on back-
%   propagation to reconstruct cavitation from selective the data collected
%   by the simulated line array using the angular spectrum method. 
%
%
% Change Log
%   201008XX   Costas Arvanitis 
%   201607XX   Calum Crake       Added frequency domain method
%   201612XX   Scott Schoen Jr   Cleanup and documentatioin; repurposing to
%                                ASA-specific processes.
%   201701XX   Scott Schoen Jr   New version to average sound speed and
%                                propagate in steps
%
% 
%**************************************************************************

% Include necessary directories
addpath( 'C:\MATLAB\inc\sjs' );

clear all;
close all;
clc;

%% Load in data file data file
disp('Loading file...');
tic;
sourceFile = ...
    '../data/results/oneBubble/layeredMedium_3000mps_04mm_1bub.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

%% Settings

% -------------------------- Plot Settings --------------------------------
% Plot layers on top of reconstruction
overPlotLayers = 0; 

% Plot sound speed profile used
plotSoundSpeedProfile = 1;

% Plot sound speed region
plotSoundSpeedRegion = 1;
% -------------------------------------------------------------------------

% --------------------- Layer Position Settings ---------------------------
% Set number of layers to divide up field into
% Start positions of layers. Set this to 0 to divide evenly into the
% desired number of layers
% layerPositions = [0, 22.5:0.5:24, 24.5:32, 32:0.5:40, 40.5:2:79 ]; % Skull 100 [mm]
% layerPositions = [25:32, 32.5:0.5:38, 39:2:81, 81.5:0.5:85.5, 86:2:103 ] - 25; % Skull 75 [mm]
% layerPositions = [0:10:30, 35:5:50, 52.5:2.5:55, 56:1:57, 58:0.4:60.5, 75]; % Stratified [mm]
layerPositions = [0, 32.1, 35.9 ]; % Layer [mm]
% layerPositions = [0, 27.1, 29.9, 34.1, 35.9 ]; % 2 Layers [mm]
numLayers = length(layerPositions);

% To use evenly spaced layers, set layerPositions = 0 and specify numLayers
% layerPositions = 0;
% numLayers = 2;
%--------------------------------------------------------------------------

% ------------------- Selective Averaging Settings ------------------------
% Determine which portion of the transverse sound speed to average 
useEntireTransverseRange = 1;
% If we want to use just portions around the target region, specify the
% point of interest and buffer (how much of transverse slice to include in
% averaging).
targetDepth = 60; % Distance from receiver array [mm]
targetTransversePosition = 0; % [mm]
transverseBuffer = 10; % [mm]
transverseSpan = 0.9; % What fraction of the array to widen to.
% -------------------------------------------------------------------------

% -------------------- Time Limit Settings --------------------------------
% Set start and stop times. If either is set to NaN, the entire range 
% is used.
t0 = 35E-6; % [s]
t1 = 85E-6; % [s]
t1 = NaN;
% -------------------------------------------------------------------------

bin = 'y';
deld = 2;
numTargets = data.numSources; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

% Medium proerties
c = 1500; % [m/s]
rho0 = 1000; % [kg/m^3]

% Set location of the probe [mm]
loc = dx*(data.yDim - 5) - data.excit_loc(numTargets,2,1)*dx;

% Get the receiver data
arrayData = data.aedata;

% Compute the time step
if ~any( isnan( [t0, t1] ) )
    t = data.t; % Time vector [s]
    
    % Get start and stop indices
    tIndex0 = find( t > t0, 1 );
    tIndex1 = find( t > t1, 1 );
    
    % Only take portion of received data between start and stop times
    t = t( tIndex0: tIndex1 );
    arrayData = arrayData( :, tIndex0 : tIndex1 );
    
end
dt = data.t(2) - data.t(1); % Time step [s]

% Get the x-positions of the sensors
y = data.yPositionVector; % [m]

% Get the size of the received US data array. This array has dimension:
%  (index of sensor) by (time index)
ss = size(arrayData);
numSensors = ss(1);
numTimeSteps = ss(2);

% Set the number of receiving channels. The number of 
channels = 64*2;

% Now resample the data for that number of sensors.
xtr = dx*(channels*floor(numSensors/channels)); % Width of "sensor" [mm]

% What are these?
rcn_x1 = 10;
rcn_x2 = numTimeSteps;

% Get the location of the left and right bounds of the field
lft_el = -xtr/2; % [mm]
rght_el = xtr/2; % [mm]
step = 0.5;%in images is 0.5;
near_coor = 0;
far_coor = dx*data.yDim;

%% Bin the data

% Get the data starting rcn_x1 steps in
rf128 = arrayData(:, rcn_x1:end );

% Initialize array to hold re-sampled data
rf = zeros( channels, numTimeSteps - rcn_x1 + 1 );

% If we're binning the data
if bin == 'y'
    
    % For each channel
    for jj = 1:channels
        
        % Find the channels of the original data that will be all grouped
        % into the new channel jj
        startChannel = (floor(numSensors/channels))*(jj-1) + 1;
        endChannel = (floor(numSensors/channels))*jj;
        
        % Assign to binned data array
        rf(jj,:) = sum( ...
            rf128(startChannel:endChannel, :) ...
            );
        
    end
    
    % Correct for different binning
    rf1 = rf./sqrt(floor(ss(1)/channels));
    
else
    
    if deld < (floor(numSensors/channels)) + 1
        for jj = 1:channels
            
            % Get sensor data to be collected together into channel jj
            startChannel = (floor(numSensors/channels))*(jj-1) + 1;
            endChannel = (floor(numSensors/channels))*(jj-1) + deld + 1;
            rf(jj,:) = sum(rf128(startChannel:endChannel,:));
        end
    else
        disp( 'deld is too large.' );
    end
    
    % Correct for different binning
    rf1 = rf./deld;
    
end

%% Get frequency parameters
centerChannelIndex = round( numSensors./2 );

Fs = 1./dt;
df = Fs./numTimeSteps;
fVector = 0 : df : ( numTimeSteps - 1 ).*df;
halfwayIndex = floor( numTimeSteps )./2;

%% Compute Field with Angular Spectrum Approach

% Set up number of padding channels (in each direction!)
padChannels = numSensors;
totalChannels = 4*128; % This is the number of chanels, plus the number of pad channels

% Set the number of points at which the field will be calculated
zEvaluationPoints = 256;

% Initialize angular spectrum map
asamap = zeros(zEvaluationPoints, totalChannels);

% Pad the array data with 0s.
paddedArrayData = padarray( double(rf), totalChannels/2 - channels/2 );

% Set the image depth
zlength = 1E-3*far_coor; % [m]

% Create array of image depths
z = linspace(0, zlength, zEvaluationPoints);
dz = z(2) - z(1);

% Create sensor position vector for padded array
channelWidth_mm = xtr*1E-3;
xlength = totalChannels*channelWidth_mm/channels;
x = linspace( -xlength/2, xlength/2, totalChannels );%for frequency sweep video

% Set up harminic frequency bins
fMin = 0.28E6; % [Hz]
fMax = 1.48E6; % [Hz]
[~,h1] = min( abs(fVector - fMin) );
[~,h2] = min( abs(fVector - fMax) );
fbins = (h1:h2)';
ss = size(fbins);

% Get time step and vector
dt = data.t(2) - data.t(1); % Time step [s]
t = data.t; % Time vector [s]

% Create frequency vector
Fs = 1./dt; % Sampling frequency [Hz]
fk = Fs.*linspace(0,1,length(rf));
tic

% Get the padded RF data into the appropriate dimension
paddedArrayData = paddedArrayData';

% Take the FFT of the padded data on each channel
paddedArrayDataTilde = fft(paddedArrayData);

% Time reconstruction
tic

% Initialize array to hold final image
asamap = zeros( length(z), totalChannels);

% Compute average sound speed at each point on the reconstruction z vector
dataC = data.c; % As a function of z *from left* (not receiver!)
dataZ = data.xPositionVector; % [m]

% Slice out the portion between the array and the greatest depth
zStartIndex = find( ... % Will fail if z depth too big
    dataZ > (data.recPosition - z(end)), 1 );
zEndIndex = find( ...
    dataZ > (data.recPosition - z(1)), 1 );

% Get the average sound speed in each plane (at each z). If the
% useEntireTransverseRange is specified, the average will be taken across
% the entirety of each slice. Otherwise, only sound speeds
if useEntireTransverseRange

    cFieldMeanPlane = mean( ...
        dataC( :, zStartIndex : zEndIndex ), 1 );
else
    
    % Initialize vector to hold averaged values
    numXPoints = length( data.yPositionVector );
    numZPoints = zEndIndex - zStartIndex + 1;
    cFieldMeanPlane = zeros( 1, numZPoints );
    
    % First define bounding lines
    xSpan = max(data.yPositionVector) - min(data.yPositionVector);
    offset = xSpan./2;
    xOffset = data.yPositionVector - offset;
    xMax = transverseSpan.*max( xOffset ); % [m]
    xMin = transverseSpan.*min( xOffset ); % [m]
    % Add buffers
    x0Top = 1E-3.*(targetTransversePosition + transverseBuffer./2); 
    x0Bot = 1E-3.*(targetTransversePosition - transverseBuffer./2);
    % Be sure we stay within data
    if x0Top > xMax
        x0Top = xMax;
    end
    if x0Bot < xMin
        x0Bot = xMin;
    end
    % Get slopes
    topSlope = ( xMax - x0Top )./(targetDepth./1E3);
    botSlope = ( xMin - x0Bot )./(targetDepth./1E3);
    % Get z vectors for before and after target z
    zSpeedRegions = fliplr( ...
        data.recPosition - dataZ( zStartIndex : zEndIndex ) );
    breakIndex = find( zSpeedRegions > targetDepth./1E3, 1 );
    z1 = zSpeedRegions( 1 : breakIndex - 1 );
    z2 = zSpeedRegions( breakIndex : end );
    % Define lines before and after target depth
    topLine = [ ...
        (x0Top - topSlope.*(z1 - targetDepth./1E3) ), ...
        (x0Top + topSlope.*(z2 - targetDepth./1E3) ) ...
        ];
    botLine = [ ...
        (x0Bot - botSlope.*(z1 - targetDepth./1E3) ), ...
        (x0Bot + botSlope.*(z2 - targetDepth./1E3) ) ...
        ];
    
    % Plot if desired
    if plotSoundSpeedRegion
        figure()
        hold all;
        [zSSP, xSSP] = meshgrid( data.xPositionVector, xOffset );
        pcolor( 1E3.*zSSP, 1E3.*xSSP, data.c );
        shading flat;

        zPlot = data.recPosition - zSpeedRegions;
        plot( 1E3.*zPlot, 1E3.*topLine, '--w' );
        plot( 1E3.*zPlot, 1E3.*botLine, '--w' );

        xlabel( 'Distance from Receiver [mm]' );
        xlim( 1E3.*[0, data.recPosition] );
        ylabel( 'Transverse Distance [mm]' );
        ylim( 1E3.*[ min(xSSP(:, 1)), max(xSSP(:, 1)) ] );
        axis equal;
        drawnow;
    end
    
    % Compute average within each region
    for zCount = 1 : numZPoints
        
        % Get limits
        xMax = topLine( zCount );
        xMin = botLine( zCount );
        % ...and the indices corresponding to those limits
        axialInd0 = find( xOffset > xMin, 1 );
        axialInd1 = find( xOffset > xMax, 1 );
        % Make sure the indices are valid
        if isempty( axialInd0 )
            axialInd0 = 1;
        end
        if isempty( axialInd1 )
            axialInd1 = numXPoints;
        end
        
        % Average over this region of the plane and store to output
        cMeanInPlane = mean( dataC( axialInd0 : axialInd1, zCount ), 1 );
        if isnan( cMeanInPlane )
           xxx = 6; 
        end
        cFieldMeanPlane(zCount) = cMeanInPlane;
    end
    
end


% Flip data, since reconstruction measures back from receiver
clippedC = fliplr( cFieldMeanPlane );

% Now just interpolate c at the values in z (instead of clippedZ )
numPointsRaw = zEndIndex - zStartIndex + 1; % Points in data
rawPoints = linspace(0, 1, numPointsRaw);
numPointsInterp = length(z); % Number of points to interpolate to
interpPoints = linspace(0, 1, numPointsInterp );
interpolatedC = interp1( rawPoints, clippedC, interpPoints );

% Initialize layer property vectors
zCenterEachLayer = zeros(1, numLayers);
cEachLayer = zeros(1, numLayers);

% Now for each frequency bin
for fCount = 1:ss(1)
    
    % Get the center frequency
    fc = fVector(fbins(fCount)); % [Hz]
    omega = 2*pi*fc; % [rad/s]
    
    % ff the data at fc
    [zz, ff] = min( abs(fk-fc) );
    
    % Get complex data at fc for each channel - CC based on single FFT outside loop
    P = paddedArrayDataTilde( ff + 1, : );
    
    % Convert P into K-space
    P = fftshift(fft(P));
    
    % Image reconstruction algorithm -------------
    dx = x(2)-x(1);
    dk = 2.*pi./dx; % Wavenumber increment
    startValue = ( -floor( length(x)/2 ) );
    endValue = ( ceil( length(x)/2 ) - 1 );
    
    % Create wavenumber vector
    k = (startValue:endValue).*dk./length(x);
    
    % Partition the sound speed field into a series of layers and perform
    % the propagation through each
    
    % Get layer positions if specified or divide up evenly
    layerPositionsSpecified = ~isequal( layerPositions, 0 );
    if layerPositionsSpecified
        numLayers = length( layerPositions ) + 1;
        zIndexLayers = zeros(1, numLayers);
        % Find index corresponding to start of each layer
        for layerCount = 1:numLayers - 1;
            zCurrentLayerStart = layerPositions( layerCount );
            layerStartIndex = find( 1E3.*z >= zCurrentLayerStart, 1 );
            zIndexLayers(layerCount) = layerStartIndex;
        end
        zIndexLayers(end) = length(z);
    else
        % Otherwise just space layers evenly
        zIndexLayers = round( (length(z)./numLayers).*( 1 : numLayers ) );
        zIndexLayers(end) = length(z); % Don't overshoot end
    end

    
    % Initialize angular spectrum at source plane
    previousAS = P;
    
    % Initialize field
    fullField = [];
    
    for layerCount = 1 : numLayers
        
        % Get start and stop indices of the axial z vector
        if layerCount == 1
            ind0 = 1;
        else
            ind0 = zIndexLayers( layerCount - 1 ) + 1;
        end
        ind1 = zIndexLayers( layerCount );
        
        % Get z vectors in that layer
        if layerCount == 1
            zLayer = z( ind0 : ind1 );
        else
            zLayer = z( ind0 : ind1 ) - z(ind0 - 1);
        end
        zvLayer = repmat(zLayer, [length(x), 1]);
        
        % Get the sound speed in that layer
        cLayer = mean( interpolatedC( ind0 : ind1 ) );
        
        % Store the center and sound speed of each layer
        zCenterIndex = floor( (ind1 + ind0)./2 );
        zCenterEachLayer( layerCount ) = z( zCenterIndex );
        cEachLayer( layerCount ) = cLayer;
        
        % Get wavenumber vectors in that layer
        kvLayer = repmat(k, [length(zLayer), 1]);
        kzLayer = sqrt( omega.^(2)./cLayer.^(2) - kvLayer.^(2) );
               
        % Now create an array of the appropriate size for the angular
        % spectrum at the previous plane
        previousASv = repmat( previousAS(end, :), [length(zLayer), 1]);
        
        % Now propagate the angular spectrum to the next layer
        nextAS = previousASv.*exp( 1j.*kzLayer.*zvLayer' );
        
        % Concatenate to the total ASA map
        fullField = [ fullField; nextAS ];
        
        % Now update "previousAS" variable
        previousAS = nextAS;
        
    end
    
    % Once we've assembled the full field, take the inverse transform to
    % get the map of the ASA
    asa = ifft( ...
            ifftshift( fullField, 2), [], 2 ...
            );
    
    % Get squared magnitude of the signal
    asa = 2*( abs(asa) ).^(2);
    
    % Add contribution of this frequency
    asamap = asamap + asa;
    
end

% Display computation time
clc;
computationTime = toc;
displayString = ...
    sprintf( 'ASA Method Computation Time: %6.2f s', computationTime );
disp( displayString );

%% Plot angular spectrum result

if plotSoundSpeedProfile
    figure()
    hold all;
    box on;
    % Interpolated sound speed profile
    layerSpeeds = plot( 1E3.*zCenterEachLayer, cEachLayer, '-or' );
    profileSpeeds = plot( 1E3.*z, interpolatedC, 'k' );
    xlabel( 'Distance from Receiver [mm]' );
    ylabel( 'Effective Sound Speed [m/s]' );
    legend( [layerSpeeds, profileSpeeds], ...
        'Layer Average', 'True Profile' );
    drawnow;
end

figure()
hold all;

% Plot reconstructed map
[ xAsaPlot, zAsaPlot ] = meshgrid( x, z );
pcolor( zAsaPlot.*1E3, xAsaPlot.*1E3, asamap );
shading flat;

% Plot layer positions
if overPlotLayers
    for layerCount = 1:numLayers
        
        % Get start and stop indices of the axial z vector
        ind0 = zIndexLayers( layerCount );
        
        % Construct plotting vectors
        leftX = 1E3.*[z(ind0), z(ind0)];
        leftY = [ -1E10, 1E10 ];
        
        % Plot
        layerLines = plot( leftX, leftY, '--', ...
            'Color', 0.6.*[1, 1, 1], ...
            'LineWidth', 1.2 );
    end
end

% Plot source positions
% Axial position relative to receiver
zSources = data.recPosition - data.sourcePositions( :, 1 );
% Lateral position relative to receiver center
xSpan = max(data.yPositionVector) - min(data.yPositionVector);
offset = xSpan./2 - dx;
xSources = data.sourcePositions( :, 2 )  - offset;
plot( zSources.*1E3, xSources.*1E3, 'ro' );

% Reverse x-direction to mirror arrangement
set( gca, 'XDir', 'reverse' );

ylim( [-40, 40] );
xlim( [0, 80] );

xlabel( 'Distance from Receiver [mm]' );
ylabel( 'Sensor Position [mm]' );

% caxis([0, 3E5]);
max(max(asamap))

cBarHandle = colorbar;

cBarHandle.Label.String = '';
cBarHandle.Label.FontSize = 16;
cBarHandle.Label.Interpreter = 'latex';
cBarHandle.TickLabelInterpreter = 'latex';

% Get axial profile
figure()
middleIndex = round(length(x)./2) + 1;
centerProfile = asamap( :, middleIndex );
centerProfileNorm = centerProfile./max(centerProfile);
plot( 1E3.*z, centerProfileNorm, 'k' );
ylim( [0, 1.01] );
xlabel( 'Distance From Receiver [mm]' );
set( gca, 'XDir', 'Reverse' );
ylabel( 'Normalized Intensity' );

% % Get axial profile
% figure()
% hold all
% middleIndex = find( x > 0, 1 );
% maxValue = max( max( asamap( :,  middleIndex-5:middleIndex + 5 ) ) );
% for cpCount = middleIndex-5:middleIndex + 5
%     centerProfile = asamap( :, cpCount )./maxValue;
%     plot( z, centerProfile + (cpCount - middleIndex).*1.5, 'k' );
% end
% set(gca, 'XDir', 'Reverse' )

