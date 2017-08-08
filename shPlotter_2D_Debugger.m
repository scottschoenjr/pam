% % 2D SPHM Debugger
% 
% % Include necessary directories
% addpath( 'C:\MATLAB\inc\sjs' );
% 
% clear all;
% close all;
% clc;
% 
% %% Load in data file data file
% disp('Loading file...');
% tic;
% sourceFile = ...
%     '../data/results/oneBubble/uniformMedium_1500mps_1bub_4040.mat';
% data = load(sourceFile);
% disp(['               ...done (', num2str(toc), ' s).' ] )
% 
% plotDelayedData = 1;
% 
% % Set spherical layer properties
% useSphericalLayer = 0;
% layerRadius = 40E-3; % [m]
% layerThickness = 4E-3; % [m]
% cLayer = 6000; % [m/s]
% rhoLayer = 1000; % [kb/m^3]
% 
% %% Bin the Recevier Data as before
% bin = 'y';
% deld = 2;
% numTargets = data.numSources; % Number of targets
% dx = data.dx*1E3; % Element spacing [mm]
% 
% % Get the receiver data
% arrayData = data.aedata;
% 
% % Compute the time step
% dt = data.t(2) - data.t(1); % Time step [s]
% 
% % Get the size of the received US data array. This array has dimension:
% %  (index of sensor) by (time index)
% ss = size(arrayData);
% numSensors = ss(1);
% numTimeSteps = ss(2);
% 
% % Set the number of receiving channels. The number of
% channels = 2*64;
% 
% % Now resample the data for that number of sensors.
% xtr = dx*(channels*floor(numSensors/channels)); % Width of "sensor" [mm]
% 
% % What are these?
% rcn_x1 = 10;
% 
% % Get the data starting rcn_x1 steps in
% rf128 = arrayData(:, rcn_x1:end );
% 
% % Initialize array to hold re-sampled data
% rf = zeros( channels, numTimeSteps - rcn_x1 + 1 );
% 
% % If we're binning the data
% if bin == 'y'
%     
%     % For each channel
%     for jj = 1:channels
%         
%         % Find the channels of the original data that will be all grouped
%         % into the new channel jj
%         startChannel = (floor(numSensors/channels))*(jj-1) + 1;
%         endChannel = (floor(numSensors/channels))*jj;
%         
%         % Assign to binned data array
%         rf(jj,:) = sum( ...
%             rf128(startChannel:endChannel, :) ...
%             );
%         
%     end
%     
%     % Correct for different binning
%     rf1 = rf./sqrt(floor(ss(1)/channels));
%     
% else
%     
%     if deld < (floor(numSensors/channels)) + 1
%         for jj = 1:channels
%             
%             % Get sensor data to be collected together into channel jj
%             startChannel = (floor(numSensors/channels))*(jj-1) + 1;
%             endChannel = (floor(numSensors/channels))*(jj-1) + deld + 1;
%             rf(jj,:) = sum(rf128(startChannel:endChannel,:));
%         end
%     else
%         disp( 'deld is too large.' );
%     end
%     
%     % Correct for different binning
%     rf1 = rf./deld;
%     
% end
% 
% %% Compute Field at the origin based on RF data
% 
% % Set the origin
% x0 = 40E-3;
% z0 = 40E-3;
% 
% % Set radius of interest (about origin)
% rVector = linspace(55E-3, 60E-3, 25); % [m]
% 
% % Get the receiver data
% arrayData = double(rf1);
% 
% % Get x (transverse) positions of sensors
% x = data.yPositionVector; % [m] (it's called y in the data creation script)
% 
% % Compute radial distance to each sensor from origin
% arrayX = linspace( min(x), max(x), channels); % Receiver positions
% dx = arrayX(2) - arrayX(1);
% arrayX = arrayX - arrayX(1); % Account for offset
% arrayZ = data.recPosition;
% 
% sensorRadialPositions = sqrt( ...
%     ( arrayZ - z0 ).^(2) + ( arrayX - x0 ).^(2) );
% 
% % Set effective radius (i.e., where the linear array signals will be
% % propagated to to simulate a circular array)
% Reff = min( sensorRadialPositions );
% 
% % Radial distance from array to target surface
% sensorRadialDifferences = sensorRadialPositions - Reff;
% 
% % Determine the minimum and maximum angles of the array
% theta = atan2( arrayX - x0, arrayZ - z0 );
% thetaMax = max(theta);
% thetaMin = min(theta);
% dTheta = (thetaMax - thetaMin)./length(theta); % Radians each element spans
% 
% % Determine the maximum m value. This is the maximum number of times a
% % signal can repeat in distance 2.*pi and still be adequately sampled
% % by dTheta spacing. Analogous to Fs = 1/dt.
% mMax = floor(2.*pi./dTheta);
% if ( mod( mMax - length(theta), 2 ) ~= 0 )
%     mMax = mMax - 1;
% end
% mVector = 0:mMax - 1;
% 
% debugDelays = 0;
% 
% if debugDelays
%     
%     % Now compute delays to propagate each sensor's data to the
%     % circle with that radius
%     testC = linspace( 1540, 1550, 100 );
%     intVariance = 0.*testC;
%     for cCount = 1:length( testC )
%         cNearArray = testC( cCount );
%         
%         % Compute delays
%         delays_seconds = sensorRadialDifferences./(cNearArray);
%         delays_samples = round(delays_seconds./dt); % No. samples to delay each channel
%         
%         % Compute amplitude corrections
%         corrections = sqrt(sensorRadialPositions./Reff);
%         
%         % Shift received data to the circle (sphere) with r = Reff
%         % TODO: There should be a way to vectorize this (shift all rows in one
%         % operation rather than looping).
%         delayedData = arrayData;
%         for channelCount = 1:channels
%             
%             % Get current channel's data
%             channelData = arrayData( channelCount, : );
%             % Get current delay
%             delaySamples = delays_samples( channelCount );
%             % Shift that channel by the negative number of samples to shift data
%             % forward
%             delayedChannelData = corrections( channelCount ).* ...
%                 circshift( channelData, [ 0, -delaySamples ] );
%             % Store back to data array
%             delayedData( channelCount, : ) = delayedChannelData;
%             
%         end
%         
%         % Chop off the end of the data, since it's just looped from the beginning
%         % of each time series (due to 'circshift')
%         maxDelaySamples = max( delays_samples );
%         lastIndex = length( data.t ) - (maxDelaySamples + 1);
%         delayedData = delayedData( :, 1:lastIndex  );
%         delayedTime = data.t( 1:lastIndex );
%         
%         
%         % Compute intensity along arraival time
%         tArrival = 46E-6;
%         tIndex = find( delayedTime > tArrival, 1 );
%         profile = delayedData( :, tIndex );
%         
%         % Store the variance of the intensity
%         intVariance( cCount ) = std( profile ).^(2);
%         
%     end
%     
%     figure()
%     plot( testC, intVariance );
%     xlabel( 'Sound Speed [m/s]' );
%     ylabel( ['Variance Along t = ', num2str( tArrival.*1E6 ), ' \mu s' ]);
% else
%     
%     % Now compute delays to propagate each sensor's data to the
%     % circle with that radius
%     testA = linspace( 0.2, 1.5, 100 );
%     intVariance = 0.*testA;
%     
%     for aCount = 1:length( testA )
%         
%         cNearArray = 1548;
%         
%         % Compute delays
%         delays_seconds = sensorRadialDifferences./(cNearArray);
%         delays_samples = round(delays_seconds./dt); % No. samples to delay each channel
%         
%         % Compute amplitude corrections
%         corrections = (sensorRadialPositions./Reff).^(testA(aCount));
%         
%         % Shift received data to the circle (sphere) with r = Reff
%         % TODO: There should be a way to vectorize this (shift all rows in one
%         % operation rather than looping).
%         delayedData = arrayData;
%         for channelCount = 1:channels
%             
%             % Get current channel's data
%             channelData = arrayData( channelCount, : );
%             % Get current delay
%             delaySamples = delays_samples( channelCount );
%             % Shift that channel by the negative number of samples to shift data
%             % forward
%             delayedChannelData = corrections( channelCount ).* ...
%                 circshift( channelData, [ 0, -delaySamples ] );
%             % Store back to data array
%             delayedData( channelCount, : ) = delayedChannelData;
%             
%         end
%         
%         % Chop off the end of the data, since it's just looped from the beginning
%         % of each time series (due to 'circshift')
%         maxDelaySamples = max( delays_samples );
%         lastIndex = length( data.t ) - (maxDelaySamples + 1);
%         delayedData = delayedData( :, 1:lastIndex  );
%         delayedTime = data.t( 1:lastIndex );
%         
%         
%         % Compute intensity along theta
%         int = sum( abs(delayedData).^(2), 2 );      
%         intVariance( aCount ) = std(int)./max( int );
%         
%     end
%     
%     figure()
%     plot( testA, intVariance );
%     xlabel( 'Amplitude Factor [m/s]' );
%     ylabel( ['Variance Along $\theta$' ]);
% end
% 
clear all
close all
clc

r0 = 0.1;
r = 0.5;

krVector = (2.*pi./1500).*linspace( 0.5E6, 15.E6, 200 ).*r;
mVector = round( linspace( 0, 50, 25 ) );
HMat = zeros( length( krVector ), length( mVector ) );

for krCount = 1:length( krVector );
    kr = krVector( krCount );
    kr0 = kr.*r0./r;
    for mCount = 1:length( mVector );
        m = mVector( mCount );
        HMat( krCount, mCount ) = besselh( m, kr )./besselh( m, kr0 );
    end
end

figure()
[krv, mv] = meshgrid( krVector, mVector );
pcolor( krv, mv, abs(HMat') );
shading flat;
xlabel( '$kr$' );
ylabel( '$m$' );
colorbar;