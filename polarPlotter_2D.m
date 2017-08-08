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
    '../data/results/oneBubble/circularLayer_6000mps_20mm_4mm_4040.mat';
% sourceFile = ...
%     '../data/analytical/data3000.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

plotDelayedData = 1;

% Set spherical layer properties
useSphericalLayer = 1;
layerRadius = 20E-3; % [m]
layerThickness = 4E-3; % [m]
cLayer = 6000; % [m/s]
rhoLayer = 1000; % [kb/m^3]

%% Bin the Recevier Data as before
bin = 'y';
deld = 2;
numTargets = data.numSources; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

% Medium proerties
c = 2000; % [m/s]
rho0 = 1000; % [kg/m^3]

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
channels = 2*64;

% Now resample the data for that number of sensors.
xtr = dx*(channels*floor(numSensors/channels)); % Width of "sensor" [mm]

% What are these?
rcn_x1 = 10;

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
Nz = 20;
Nx = 20;
zGuess = linspace( 35./1E3, 65./1E3, Nz ); % [m]
xGuess = linspace( 35./1E3, 45./1E3, Nx ); % [m]
pMap = zeros( Nz, Nx );

tic;
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
        dx = arrayX(2) - arrayX(1);
        arrayX = arrayX - dx; % Account for offset
        arrayZ = data.recPosition;
        
        sensorRadialPositions = sqrt( ...
            ( arrayZ - origin(1) ).^(2) + ( arrayX - origin(2) ).^(2) );
        
        % Radial distance from array to target surface
        sensorRadialDifferences = sensorRadialPositions - originToArray;
        
        %%%%%%% DEBUG %%%%%%%
        testR = mean( sensorRadialPositions );
        r1 = testR;
        
        sensorRadialPositions = sqrt( ...
            ( arrayZ + dx./2 - origin(1) ).^(2) + ( arrayX + dx./2 - origin(2) ).^(2) );
        
        % Radial distance from array to target surface
        sensorRadialDifferences = sensorRadialPositions - testR;
        %%%%%%%%%%%%%%%%%%%%%
        
        
        % Define angle
        theta = acos( r1./sensorRadialPositions );
        
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
        
        % Plot if desired
        if plotDelayedData && xCount == 15 && zCount == 15
            
            % Plot delayed data
            figure();
            subplot( 2, 1, 1 );
            pcolor( delayedData );
            shading flat;
            subplot( 2, 1, 2 );
            pcolor( arrayData );
            shading flat;
            
            % Plot energy as a function of angle
            figure()
            hold all;
            theta = atan2( arrayX - origin(2), arrayZ - origin(1) );
            energy = sum( delayedData.^(2), 2 );
            plot( 180.*theta./(2.*pi), energy./max(abs(energy)) );
            ylabel( 'Normalized Energy' );
            xlabel( 'Angle [deg]' );
            drawnow;
            
        end
        
        % Chop off the end of the data, since it's just looped from the beginning
        % of each time series (due to 'circshift')
        maxDelaySamples = max( delays_samples );
        lastIndex = length( data.t ) - (maxDelaySamples + 1);
        %         delayedData = delayedData( :, 1:lastIndex  );
        %         delayedTime = data.t( 1:lastIndex );
        delayedData( :, lastIndex:end  ) = 0;
        delayedTime = data.t;
        
        % Now, take the FFT of the data. Zero pad the data in time to the next
        % power of 2 (where FFT is most efficient).
        numTimeSteps = length( delayedTime );
        numSamplesWithPadding = 2.^( nextpow2(numTimeSteps) ); % <--- CHECK THIS!
        delayedDataTilde = fft( delayedData, numTimeSteps, 2 ); % Along second dim.
        
        % Get frequency parameters
        Fs = 1./dt;
        df = Fs./numTimeSteps;
        fVector = 0 : df : ( numTimeSteps - 1 ).*df;
        halfwayIndex = floor( numTimeSteps )./2;
        
        % Set bandwidth of receiver
        fMin = 0.28E6;
        fMax = 1.5E6;
        
        % Compute transmission coefficient for shell
        
        %% Loop over frequencies to propagate the signal to the origin
        pAtOriginTilde = 0.*delayedDataTilde;
        r0 = dx./sqrt(2); % Reconstruction radius is pixel radius
        for fCount = 1:length( fVector )
            
            % Get the current frequency
            f = fVector(fCount);
            
            % Make sure we're in the frequency range of interest
            fInRange = ( f >= fMin ) && ( f <= fMax );
            if ~fInRange
                continue;
            end
            
            % Compute the transfer function
            omega = 2.*pi.*f;
            k = omega./c1;
            
            if useSphericalLayer
                
                % Get distance to outside of layer
                rIn = layerRadius;
                rOut = layerRadius + layerThickness;
                
                % Note that omega must be a scalar, so this computation has
                % to be within the loop
                T = ...
                    shellTransmission( ...
                    1500, 1000, ... % c, rho Inside shell
                    3000, 1000, ... % c, rho Within shell
                    1500, 1000, ... % c, rho Outside shell
                    rIn, rOut, omega );
                
                % Propagation from synthetic "curved" array to layer
                H1 = (r1./rOut).*exp( -1i.*k.*(r1 - rOut) );
                
                % Propagate through layer
                H2 = 1./T;
                
                % Propagate from layer to test point
                H3 = (rIn./r0).*exp( -1i.*k.*(rIn - r0) );
                
                % Compute total transfer function
                H = H3.*H2.*H1;
                
                %%%%%% DEBUG %%%%%%%
                if fCount == 255
                    TMat( zCount, xCount ) = T;
                end
                Hmat( zCount, xCount ) = H;
                %%%%%%%%%%%%%%%%%%%
                
            else
                H = (r1./r0).*exp( -1i.*k.*(r1 - r0) );
                
                %%%%%% DEBUG %%%%%%%
                Hmat( zCount, xCount ) = H;
                %%%%%%%%%%%%%%%%%%%
                
            end
            
            % Back propagate the signal to the assumed origin
            pTilde = delayedDataTilde(:, fCount);
            pAtOriginTilde(:, fCount) = H.*pTilde;
            
        end
        
        % Get the intensity contribution of th
        term1 = abs( sum( pAtOriginTilde ) ).^(2);
        term2 = sum( abs(delayedDataTilde).^(2) );
        
        % Get the total power from that pixel (integrate over time)
        pressureAmplitude = sum( term1 - term2 );
        
        % Assign intensity to that pixel
        pMap( zCount, xCount ) = pressureAmplitude;
        
    end
end
reconTime = toc;
disp([ 'Reonstruction Time: ', num2str(reconTime) ]);

%% Plot reconstructed image
figure()
hold all;

[xPlot, zPlot ] = meshgrid( xGuess, zGuess );
pcolor( zPlot.*1E3, xPlot.*1E3, abs( pMap ) );
plot( 1E3.*data.sourcePositions(1), 1E3.*data.sourcePositions(2), 'ro' );
shading flat;
xlabel( 'Axial Distance [mm]' );
ylabel( 'Transverse Distance [mm]' );

% Plot magnitude and phase of transfer function
figure()
pcolor( zPlot.*1E3, xPlot.*1E3, abs( Hmat ) );
shading flat;
xlabel( 'Axial Distance [mm]' );
ylabel( 'Transverse Distance [mm]' );

figure()
pcolor( zPlot.*1E3, xPlot.*1E3, angle( Hmat ) );
shading flat;
xlabel( 'Axial Distance [mm]' );
ylabel( 'Transverse Distance [mm]' );

% Plot axial profile
 centerIndex = round( length( pMap(:, 2) )./2 );
profile = abs(pMap(:, centerIndex) );
fwhm( profile./max(profile), zGuess.*1E3, 1 );
xlabel( 'Axial Distance [mm]' );
ylabel( 'Normalized Intensity' );



