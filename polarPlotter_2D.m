%**************************************************************************
%
% Polar Field Reconstruction (2D)
%
%   This script uses a passive acoustic mapping (PAM) based on back-
%   propagation of spherical (circular) waves. 
%
%
% Change Log
%   20170216   Scott Schoen Jr   Initial Version
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
    '../data/results/oneBubble/uniformMedium_1500mps_1bub.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )



%% Bin the Recevier Data as before
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
dt = data.t(2) - data.t(1); % Time step [s]

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

%% Compute Field at the origin based on RF data

% Get the receiver data
arrayData = double(rf);

% Get x (transverse) positions of sensors
x = data.yPositionVector; % [m] (it's called y in the data creation script)

% Set the range of interrogation points
Nz = 50;
Nx = 50;
zGuess = linspace( 38./1E3, 42./1E3, Nz ); % [m]
xGuess = linspace( 38./1E3, 42./1E3, Nx ); % [m]
pMap = zeros( Nz, Nx );

for xCount = 1:Nx
    for zCount = 1:Nz
        
        % Now set the origin point (assumed location of the source)
        origin = [ zGuess(zCount), xGuess(xCount) ]; % (z, x) [m]
        
        % Get the z-distance from the origin to the nearest array element
        originToArray = abs( data.recPosition - origin(1) );
        r1 = originToArray;
        
        % Now compute delays to propagate each sensor's data to the 
        % circle with that radius
        
        % Get sound speed near array
        % TODO: read this from the data
        cNearArray = 1500;
        c1 = cNearArray;
        
        % Compute radial distance to each sensor from origin
        arrayX = linspace( min(x), max(x), channels); % Receiver positions
        arrayZ = data.recPosition;
        
        sensorRadialPositions = sqrt( ...
            ( arrayZ - origin(1) ).^(2) + ( arrayX - origin(2) ).^(2) );
        
        % Radial distance from array to target surface
        sensorRadialDifferences = sensorRadialPositions - originToArray;
        
        % Compute delays
        delays_seconds = sensorRadialDifferences./cNearArray;
        delays_samples = round(delays_seconds./dt); % No. samples to delay each channel
        
        % Shift received data to the circle (sphere) that just touches the array.
        % TODO: There should be a way to vectorize this (shift all rows in one
        % operation rather than looping).
        delayedData = arrayData;
        for channelCount = 1:channels
            
            % Get current channel's data
            channelData = arrayData( channelCount, : );
            % Get current delay
            delaySamples = delays_samples( channelCount );
            % Shift that channel by the negative number of samples to shift data
            % forward
            delayedChannelData = circshift( channelData, [ 0, -delaySamples ] );
            % Store back to data array
            delayedData( channelCount, : ) = delayedChannelData;
            
        end
        
        % Chop off the end of the data, since it's just looped from the beginning
        % of each time series (due to 'circshift')
        maxDelaySamples = max( delays_samples );
        lastIndex = length( data.t ) - (maxDelaySamples + 1);
        delayedData = delayedData( :, 1:lastIndex  );
        delayedTime = data.t( 1:lastIndex );
        
        % Now, take the FFT of the data. Zero pad the data in time to the next
        % power of 2 (where FFT is most efficient).
        numTimeSteps = length( delayedTime );
        numSamplesWithPadding = 2.^( nextpow2(numTimeSteps) );
        delayedDataTilde = fft( delayedData, numTimeSteps, 2 ); % Along second dim.
        
        % Get frequency parameters
        Fs = 1./dt;
        df = Fs./numTimeSteps;
        fVector = 0 : df : ( numTimeSteps - 1 ).*df;
        halfwayIndex = floor( numTimeSteps )./2;
        
        % Set bandwidth of receiver
        fMin = 0.28E6;
        fMax = 1.5E6;
        
        %% Loop over frequencies to propagate the signal to the origin
        pAtOrigin = 0;
        for fCount = 1:length( fVector )
            
            % Get the current frequency
            f = fVector(fCount);
            
            % Make sure we're in the frequency range of interest
            fInRange = ( f >= fMin ) && ( f <= fMax );
            if ~fInRange
                continue;
            end
            
            % Get the propagation factor
            omega = 2.*pi.*f;
            k = omega./c1;
            propagationFactor = r1.*exp(-1i.*omega.*r1);
            
            % Back propagate the signal to the assumed origin
            pTilde = delayedDataTilde(:, fCount);
            pAtOrigin = pAtOrigin + pTilde.*propagationFactor;
            
        end
        
        pMap( zCount, xCount ) = sum( pAtOrigin.^(2) );
        
    end
end

% Plot
figure()

[zPlot, xPlot] = meshgrid( zGuess, xGuess );
pcolor( zPlot.*1E3, xPlot.*1E3, abs( pMap ) );
shading flat;
xlabel( 'Axial Distance [mm]' );
ylabel( 'Transverse Distance [mm]' );

