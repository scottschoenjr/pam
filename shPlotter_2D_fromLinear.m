%**************************************************************************
%
% Circular Field Reconstruction (2D)
%
%   This script uses a passive acoustic mapping (PAM) based on back-
%   propagation of spherical (circular) waves. This is a 2D, radially
%   symmetric version of the spherical harmonic transform.
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
    '../data/results/oneBubble/uniformMedium_1500mps_1bub_4040_circularArray_45.mat';
% sourceFile = ...
%     '../data/results/threeBubble/uniformMedium_1500mps_3bub.mat';
% sourceFile = ...
%     '../data/analytical/data3000.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

plotDelayedData = 1;

% Set spherical layer properties
useSphericalLayer = 0;
layerRadius = 40E-3; % [m]
layerThickness = 4E-3; % [m]
cLayer = 6000; % [m/s]
rhoLayer = 1000; % [kb/m^3]

%% Bin the Recevier Data as before
bin = 'y';
deld = 2;
numTargets = data.numSources; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

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

% Set the origin
x0 = 40E-3;
z0 = 40E-3;

% Set radius of interest (about origin)
rVector = linspace(10E-3, 60E-3, 25); % [m]

% Get the receiver data
arrayData = double(rf1);

% Get x (transverse) positions of sensors
x = data.yPositionVector; % [m] (it's called y in the data creation script)

% Compute radial distance to each sensor from origin
arrayX = linspace( min(x), max(x), channels); % Receiver positions
dx = arrayX(2) - arrayX(1);
arrayX = arrayX - arrayX(1); % Account for offset
arrayZ = data.recPosition;

sensorRadialPositions = sqrt( ...
    ( arrayZ - z0 ).^(2) + ( arrayX - x0 ).^(2) );

% Set effective radius (i.e., where the linear array signals will be
% propagated to to simulate a circular array)
Reff = min( sensorRadialPositions );

% Radial distance from array to target surface
sensorRadialDifferences = sensorRadialPositions - Reff;

% Determine the minimum and maximum angles of the array
theta = atan2( arrayX - x0, arrayZ - z0 );
thetaMax = max(theta);
thetaMin = min(theta);
dTheta = (thetaMax - thetaMin)./length(theta); % Radians each element spans

% Determine the maximum m value. This is the maximum number of times a
% signal can repeat in distance 2.*pi and still be adequately sampled
% by dTheta spacing. Analogous to Fs = 1/dt.
mMax = floor(2.*pi./dTheta);
if ( mod( mMax - length(theta), 2 ) ~= 0 )
    mMax = mMax - 1;
end
mVector = 0:mMax - 1;

tic;

% % Now compute delays to propagate each sensor's data to the
% % circle with that radius
%
% % Get sound speed near array
% % TODO: read this from the data
% cNearArray = 1548;
% c1 = cNearArray;
%
% % Compute delays
% delays_seconds = sensorRadialDifferences./(cNearArray);
% delays_samples = round(delays_seconds./dt); % No. samples to delay each channel
%
% % Compute amplitude corrections
% corrections = sqrt(sensorRadialPositions./Reff);
% corrections = (sensorRadialPositions./Reff).^(0.8);
%
% % Shift received data to the circle (sphere) with r = Reff
% % TODO: There should be a way to vectorize this (shift all rows in one
% % operation rather than looping).
% delayedData = arrayData;
% for channelCount = 1:channels
%
%     % Get current channel's data
%     channelData = arrayData( channelCount, : );
%     % Get current delay
%     delaySamples = delays_samples( channelCount );
%     % Shift that channel by the negative number of samples to shift data
%     % forward
%     delayedChannelData = corrections( channelCount ).* ...
%         circshift( channelData, [ 0, -delaySamples ] );
%     % Store back to data array
%     delayedData( channelCount, : ) = delayedChannelData;
%
% end
%
% % Chop off the end of the data, since it's just looped from the beginning
% % of each time series (due to 'circshift')
% maxDelaySamples = max( delays_samples );
% lastIndex = length( data.t ) - (maxDelaySamples + 1);
% delayedData = delayedData( :, 1:lastIndex  );
% delayedTime = data.t( 1:lastIndex );

%%%%%%% DEBUG %%%%%%%
delayedData = arrayData;
delayedTime = data.t;
c1 = 1500;
%%%%%%%%%%%%%%%%%%%%%

% Plot delayed data
plotDelayedData = 0;
if plotDelayedData
    
    figure()
    int = sum( abs(delayedData).^(2), 2 ); % - sum( abs(delayedData).^(2), 2 );
    plot( 180.*theta./pi, int );
    xlabel( '$\theta$ [deg]' );
    ylabel( 'Intensity' );
    
    figure()
    [t, n] = meshgrid( delayedTime, 1:channels );
    pcolor( t.*1E6, n, delayedData );
    shading flat;
    xlabel( 'Time [ms]' );
    xlim( [20, 50] );
    ylabel( 'Channel Number' );
    return;
end

% Pad data in space
padLength = floor( (mMax - length(theta))./2 );
delayedData = padarray( delayedData, [padLength, 0] );

% Now, take the FFT of the data (in time)
numTimeSteps = length( delayedTime );
delayedDataTilde = fft( delayedData, numTimeSteps, 2 ); % Along second dim.

% Get the spherical spectrum (really circular since 2D case)
S0 = fft( delayedDataTilde );

% Get frequency parameters
Fs = 1./dt;
df = Fs./numTimeSteps;
fVector = 0 : df : ( numTimeSteps - 1 ).*df;
halfwayIndex = floor( numTimeSteps )./2;

%%%%%%% DEBUG %%%%%%%
[mPlot, fPlot] = meshgrid( mVector, fVector );
pcolor( fPlot./1E6, mPlot, abs(S0)'./max(max(abs(S0))));
xlabel('Frequency [MHz]');
xlim( [0, Fs./2]./1E6 );
ylim( [0, mMax./2] );
ylabel('$m$');
shading interp;
colorbar;
%%%%%%%%%%%%%%%%%%%%

% Set bandwidth of receiver
fMin = 0.28E6;
fMax = 1.5E6;

%%%%%% DEBUG %%%%%%%
counter = 1;
%%%%%%%%%%%%%%%%%%%%

% Compute intensity along theta at each radius
sMap = zeros( length( rVector ), mMax );
for rCount = 1:length( rVector );
    
    rFinal = rVector( rCount );
    
    % Loop over frequencies to propagate the signal to desired radius
    pReconTilde = 0.*delayedDataTilde;
    Hvec = 0.*fVector;
    
    for fCount = 1:length( fVector )
        
        % Get the current frequency
        f = fVector(fCount);
        
        % Get the spectrum at that frequency
        currentS0 = S0( :, fCount );
                
        % Make sure we're in the frequency range of interest
        fInRange = ( f >= fMin ) && ( f <= fMax );
        if ~fInRange
            continue;
        end
        
        %%%%%%%% DEBUG %%%%%%%
        if rCount == 1
            figure(1)
            currentP0 = abs(delayedDataTilde( :, fCount ));
            th = dTheta:dTheta:2.*pi - dTheta;
            subplot( 2, 1, 1)
            plot( mVector, abs(currentS0)./max(max( abs(S0) )) );
            xlabel( '$m$' );
            xlim( [0, mMax./2] );
            ylabel( '$|\mathcal{S}|$' );
            ylim([0, 1.05]);
            title( ['$f =$ ',num2str( f./1E6 ), ' MHz' ] );
            subplot( 2, 1, 2 );
            plot( th, currentP0 )
            set( gca, ...
                'XTick', [0, pi./2, pi, 3.*pi./2, 2.*pi], ...
                'XTickLabel', ...
                    {'$0$';'$\pi/2$';'$\pi$';'$3\pi/2$';'$2\pi$'} );
            xlabel( '$\theta$' );
            xlim( [0, 2.*pi] );
            ylim( [0, 3] );
            ylabel('$\tilde{p}$');
            drawnow;
            M( counter ) = getframe();
            counter = counter + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % Compute the transfer function
        omega = 2.*pi.*f;
        k = omega./c1;
        
        % Compute transfer function for that radius (function of m)
        Hnum = besselh( mVector, k.*rFinal );
        Hden = besselh( mVector, k.*Reff );
        H = (Hnum./Hden)';
        
        % Retain only the first few ms if we are back-propagating. Because
        % H(kr) dies off quickly with kr for higher m, keeping these terms
        % would cause likely noise to be amplified substantially.
        if rFinal < Reff
            H( 15:end ) = 0;
        end
        
        % Apply transfer function to spherical spectrum
        S = H.*currentS0;
        
        % Take IFFT to get presure (still in freq. domain) at new radius
        pressureAtR = ifft( S );
        pReconTilde( :, fCount ) = pressureAtR;
        
    end
    
    % Get time domain signal
    pRecon = ifft( pReconTilde, [], 2 );
    
    % Get the total energy over time for that sensor
    intensity = sum( abs( real( pRecon ) ), 2 );
    
    % Save the intensity at that radius
    sMap( rCount, : ) = intensity;
    
end

reconTime = toc;
disp([ 'Reonstruction Time: ', num2str(reconTime) ]);

%% Plot pressure at that radius
figure()
hold all;

startIndex = padLength + 1;
endIndex = length( intensity ) - padLength;
plot( theta.*180./pi, intensity(startIndex:endIndex)./max(abs(intensity(startIndex:endIndex))) );

shading flat;
xlabel( 'Azimuth [deg]' );
ylabel( 'Intensity [AU]' );
ylim( [0.5, 1] );

% Retain only the real channels of S
startIndex = padLength + 1;
endIndex = length( S0( :, 1 ) ) - padLength;
sMap = sMap( :, startIndex : endIndex );

figure()
hold on;

[rPlot, thetaPlot] = meshgrid( rVector, theta );
zPlot = rPlot.*cos( thetaPlot );
xPlot = rPlot.*sin( thetaPlot );
pcolor( zPlot, xPlot, sMap' );
shading flat;

% Get gradient on re-sampled grid

% Resample to make plot less busy
factor = 5; % Factor to reduce num points by
rStep = factor;
thetaStep = factor;
% Create vectors of indices to keep
[rDim, thetaDim] = size( sMap );
rKeepIndices = rStep : rStep : rDim - rStep;
thetaKeepIndices = thetaStep : thetaStep : thetaDim - thetaStep;
% Create interpolation points
[ rInterp, thetaInterp ] = meshgrid( ...
    rVector( rKeepIndices ), theta( thetaKeepIndices ) );
zInterp = rInterp.*cos(thetaInterp);
xInterp = rInterp.*sin(thetaInterp);
% Interpolate function
sMapResampled = interp2( rPlot, thetaPlot, sMap', ...
    rInterp, thetaInterp );

% Get the gradient
[zComp, xComp] = gradient( sMapResampled );

% Plot!
quiver( zInterp, xInterp, zComp, xComp, ...
    'Color', 'w');

axis equal;
xlabel( 'Axial Distance [mm]' );
ylabel( 'Transverse Distance [mm]' );
