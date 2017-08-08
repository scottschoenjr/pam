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
    '../data/results/oneBubble/uniform/uniformMedium_1500mps_1bub_4040_circularArray.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

% Plotting options --------------------------------------------------------
plotDelayedData = 0; % (1 to plot, 2 to plot and return)
plotSphericalSpectrum = 0;
shiftS0 = 1;
% To plot the spherical spectrum (over the frequency range) at a given
% radius
plotSpectrumAtRadius = 0;
rOfInterest = 45E-3;
% To plot the transfer function and reconstructed pressure at each radius
plotTransferFunction = 0;
plotEvery = 10;
% Plot reconstructed pressure at each radius
plotPRecon = 1;
% Plot gradient on top of reconstructed field
plotGradient = 0;
%-------------------------------------------------------------------------

% Set radius of interest (about origin)
rVector = linspace(55E-3, 22E-3, 12); % [m]

% Set the number of receiving channels.
channels = 128;

% Set the padding length
padLength = 60;

% Set number of spatial points to keep
maxMToKeep = 5;

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

% Get the receiver data
arrayData = double(rf1);

% Rename it to "delayed data" since we could get the pressure on a circular
% array by delaying the results from a linear array.
delayedData = arrayData;
delayedTime = data.t;
c1 = 1500;

% Cut off signal early to eliminate high frequency variability in the
% spectrum
tMin = 0E-6;
tMax = 60E-6;
tMaxIndex = find ( delayedTime > tMax, 1 );
tMinIndex = find ( delayedTime > tMin, 1 );
delayedData = delayedData( :, tMinIndex : tMaxIndex );
delayedTime = delayedTime( tMinIndex : tMaxIndex );

% Get x (transverse) positions of sensors
x = data.yPositionVector; % [m] (it's called y in the data creation script)

% Get the radius on which the signal was reocorded
Reff = data.R0;

% Get the angular span of the receiver array
ySpan = max( data.yPositionVector ) - min( data.yPositionVector );
thetaMax = asin( (ySpan./2)./Reff );
thetaMin = asin( -(ySpan./2)./Reff );

% Create vectors
theta = linspace( thetaMin, thetaMax, channels );
dTheta = abs(theta(2) - theta(1)); % [rad]

% Determine the maximum m value. This is the maximum number of times a
% signal can repeat in distance 2.*pi and still be adequately sampled
% by dTheta spacing. Analogous to Fs = 1/dt.
mMax = floor(2.*pi./dTheta);
if ( mod( mMax - length(theta), 2 ) ~= 0 )
    mMax = mMax - 1;
end
mVector = ( (-mMax/2) : (mMax/2) - 1 ); % m is even, so this is okay
% mVector = 0:mMax - 1;

tic;

%%%%%%% DEBUG %%%%%%%%%
% Extend theta to the extra (virtual) channels
newThetaMin = thetaMin - padLength.*dTheta;
newThetaMax = thetaMax + padLength.*dTheta;
theta = newThetaMin : dTheta : newThetaMax;
% Make sure theta has even number of points
if mod( length( theta ), 2 ) ~= 0
    theta = [ newThetaMin - dTheta, theta ];
end
% Assemble new m vector
mMax = length( theta );
mVector = -mMax./2 : mMax./2 - 1;
disp( [...
    'PADDING: Padding array data with ', ...
    num2str(padLength), ' channels on each side (symmetric).'] );
%%%%%%%%%%%%%%%%%%%%%%%
delayedData = padarray( delayedData, [padLength, 0], 'symmetric' );

%%%%%% DEBUG %%%%%%%%%%%
% Apply apodization function
wN = length( theta );
alpha = 2.*(1 + 0.0).*padLength./wN; % Fraction of channels apodized
[~, tN] = size( delayedData ); % Number of time steps
windowFunction = windowFunc( wN, 'Tukey', alpha );
windowFunction = repmat( windowFunction, tN, 1)'; % To allow matrix multiplication
delayedData = windowFunction.*delayedData;
disp( [...
    'WINDOW: Applying Tukey window to padded data'] );
%%%%%%%%%%%%%%%%%%%%%%%%

% Take the FFT of the data (in time)
numTimeSteps = length( delayedTime );
delayedDataTilde = fft( delayedData, numTimeSteps, 2 ); % Along second dim.

% Get frequency parameters
Fs = 1./dt;
df = Fs./numTimeSteps;
fVector = 0 : df : ( numTimeSteps - 1 ).*df;
halfwayIndex = floor( numTimeSteps )./2;

%%%%%%% DEBUG %%%%%%%
% Plot delayed data
if plotDelayedData
    
    figure(1000)
    power = sum( delayedData.^(2), 2 );
    [t, n] = meshgrid( delayedTime, 1:channels + 2.*padLength );
    pcolor( t.*1E6, n, delayedData );
    shading flat;
    xlabel( 'Time [ms]' );
    ylabel( 'Channel Number' );
    title( 'Received Data' );
    
    figure(1001)
    plot( linspace( 0, 360, length(mVector)), power./max(power), 'k' );
    xlabel('$\theta$ [deg]' );
    ylabel('Normalized Power' );
    title( 'Power at Receiver' );
    
    figure(1002)
    [fPlot, cPlot] = meshgrid( fVector, 1:channels+2*padLength );
    pcolor( fPlot./1E6, cPlot, abs(delayedDataTilde) );
    xlabel('Frequency [MHz]');
    xlim( [0, Fs./2]./1E6 );
    ylim( [0, mMax] );
    ylabel('Channel No.');
    title( 'Received Data Spectrum' );
    shading flat;
    cbh = colorbar;
    ylabel( cbh, '$|\tilde{p}(\omega)|$', ...
        'FontSize', 24, ...
        'interpreter', 'LaTeX' );
    
    if plotDelayedData == 2
        return;
    end
end
%%%%%%%%%%%%%%%%%

% % Filter
% fCutoff = 0.5;
% [bn, an] = butter( 15, fCutoff, 'low' );
% delayedDataTilde = filter( bn, an, delayedDataTilde, [], 2 );

% % Filter the data along theta. Numerical noise in the recorded data leads
% % to problems in getting the spherical spectrum. This low-pass filter will
% % leave only larger-scale variations. (This is essentially to filter noise)
% mCutoff = 3; % Max number of cycles between thetaMin and thetaMax
% Ft = channels./(thetaMax - thetaMin); % Samples per radian
% [bn, an] = butter( 15, mCutoff./(Ft./2) );
% delayedDataTilde = filter( bn, an, delayedDataTilde, [], 2 );

% Get the spherical spectrum (really circular since 2D case) and shift to
% align with center channel
if shiftS0
    S0 = fftshift( fft( delayedDataTilde ), 1 );
else
    S0 = fft( delayedDataTilde, [], 1 );
end

%%%%%%% DEBUG %%%%%%%
if plotSphericalSpectrum
    figure(1004)
    [mPlot, fPlot] = meshgrid( mVector, fVector );
    pcolor( fPlot./1E6, mPlot, abs(S0)'./max(abs(S0(:))) );
    xlabel('Frequency [MHz]');
    %     xlim( [0, Fs./2]./1E6 );
    %     ylim( [-mMax./2, mMax./2] );
    ylabel('$m$');
    shading interp;
    cbh = colorbar;
    ylabel( cbh, '$|\mathcal{S}_{\tilde{p},0}(m)|$', ...
        'FontSize', 24, ...
        'interpreter', 'LaTeX' );
    if plotSphericalSpectrum == 2
        return;
    end
end
%%%%%%%%%%%%%%%%%%%%

% Set bandwidth of receiver
fMin = 0.28E6;
fMax = 1.5E6;

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
        
        % Get the spherical spectrum at that frequency (funtion of m)
        currentS0 = S0( :, fCount );
        
        % Make sure we're in the frequency range of interest
        fInRange = ( f >= fMin ) && ( f <= fMax );
        if ~fInRange
            continue;
        end
        
        %%%%%%%% DEBUG %%%%%%%
        if plotSpectrumAtRadius
            
            % Evaluate only at the radius of interst
            [~, rIndex] = min( abs( rVector - rOfInterest ) );
            if rCount ~= rIndex
                continue;
            end
            
            find( rVector > rOfInterest, 1 );
            figure(998)
            currentP0 = abs(delayedDataTilde( :, fCount ));
            th = theta + pi;
            
            subplot( 3, 1, 1 );
            plot( th, currentP0./max(abs(delayedDataTilde(:))), 'k' );
            set( gca, ...
                'XTick', [0, pi./2, pi, 3.*pi./2, 2.*pi], ...
                'XTickLabel', ...
                {'$0$';'$\pi/2$';'$\pi$';'$3\pi/2$';'$2\pi$'} );
            xlabel( '$\theta$ [rad]', 'FontSize', 24 );
            xlim( [0, 2.*pi] );
            ylim( [0, 1] );
            ylabel('$\tilde{p}(r, \theta)$', 'FontSize', 24);
            title( ['$f =$ ',num2str( f./1E6 ), ' MHz' ] );
            
            subplot( 3, 1, 2)
            plot( mVector, abs(currentS0)./max( abs( S0(:)) ), 'k' );
            xlim( [0, mMax./2] );
            set( gca, 'XTickLabel', '' );
            ylabel( '$|\mathcal{S}_{\tilde{p}}(r, m)|$', 'FontSize', 24 );
            ylim([0, 1.05]);
            
            subplot( 3, 1, 3)
            plot( mVector, angle(currentS0), 'k' );
            xlabel( '$m$', 'FontSize', 24 );
            xlim( [0, mMax./2] );
            ylabel( '$\angle\,\mathcal{S}_{\tilde{p}}(r, m)$', 'FontSize', 24 );
            ylim([-pi, pi]);
            
            drawnow;
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % Compute the transfer function
        omega = 2.*pi.*f;
        k = omega./c1;
        
        % Compute transfer function for that radius (function of m)
        Hnum = besselh( mVector, k.*rFinal );
        Hden = besselh( mVector, k.*Reff );
        H = (Hnum./Hden)';
        
        % We'll need to filter in k-space, per Sec. 5.2.2 of Ref. [1].
        % According to Eq. (5.19) of that text, the max m value is
        %    |m| < pi*r/( 23.3 d/D ),
        % where d = r0 - r and D is the dynamic range. Setting D = 100 (dB)
        % gives
        D = 5E-5;
        a = min(rVector);
        d = Reff - min(rVector);
        if ~exist( 'maxMToKeep', 'var' ) || ...
                isnan( maxMToKeep ) || isequal( maxMToKeep, 0 )
            maxMToKeep = ceil( pi.*a./( 23.3.*D./d ) );
            disp( [...
                'Using ', ...
                num2str(2.*maxMToKeep + 1),...
                ' Spatial Points (m).' ...
                ] );
        end
        mKeepStart = -maxMToKeep;
        mKeepEnd = maxMToKeep;
        mStartIndex = find( mVector > mKeepStart, 1) - 1;
        mEndIndex = find( mVector > mKeepEnd, 1) - 1;
        keepIndices = zeros( 1, length(H) );
        keepIndices( mStartIndex : mEndIndex ) = 1;
        discardIndices = find( ~keepIndices );
        if rFinal < Reff
            H( discardIndices ) = 0;
        end
        
        % Apply transfer function to spherical spectrum
        S = H.*currentS0;
        
        % Take IFFT to get presure (in freq. domain) at new radius
        if shiftS0
            pressureAtR = ifft( ifftshift( S, 1 ) );
        else
            pressureAtR = ifft( S, [], 1 );
        end
        pReconTilde( :, fCount ) = pressureAtR;
        
        %%%%%%% DEBUG %%%%%%%
        if plotTransferFunction && mod( fCount, plotEvery ) == 0;
            figure(997)
            plot( mVector, abs(H) );
            ylim( [0, 4] );
            xlabel( '$m$' );
            ylabel( '$|H(m)|$' );
            title( ['$r$ = ', num2str(rVector( rCount ).*1E3 ), ' mm, '...
                '$f$ = ', num2str(f./1E6 ), 'MHz'] );
            
            figure(996)
            hold all;
            % Delete previous lines
            delete( findobj( 'Tag', 'OldLine' ) );
            normValue = max(2.*pi.*abs(delayedDataTilde(:)));
            % Plot data that's being manipulated
            plot( theta  + pi, ...
                abs( pressureAtR )./normValue, '--k', ...
                'Tag', 'OldLine' );
            % Plot over the top the actual data
            keepInds = find( theta > thetaMin & theta < thetaMax );
            plot( theta(keepInds) + pi, ...
                abs( pressureAtR(keepInds) )./normValue, 'k', ...
                'Tag', 'OldLine' );
            ylim( [0, 1] );
            title( ['$r$ = ', num2str(rVector( rCount ).*1E3 ), ' mm, '...
                '$f$ = ', num2str(f./1E6 ), 'MHz'] );
            xlabel( '$\theta$' );
            ylabel( 'Normalized $|\tilde{p}(\theta)|$' );
            
            drawnow;
        end
        %%%%%%%%%%%%%%%%%%%%%
        
    end
    
    % Get time domain signal
    pRecon = ifft( pReconTilde, [], 2 );
    
    %%%%%%% DEBUG %%%%%%%
    % Plot reconstructed pressure (i.e., the time-series pressure that
    % would be measured at that radius)
    if plotPRecon;
        figure(995)
        [tp, cp] = meshgrid( delayedTime, 1:channels + 2.*padLength );
        pcolor( tp.*1E6, cp, real( pRecon ) );
        shading flat;
        xlabel( 'Time [$\mu$s]' );
        ylabel( 'Channel Number' );
        title( ['$r$ = ', num2str(rVector( rCount ).*1E3 ), ' mm'] );
        caxis( [-30, 30] );
        drawnow;
        pause(0.5);
    end
    %%%%%%%%%%%%%%%%%%%%%
    
    % Get the total energy over time for that sensor
    intensity = sum( abs( real( pRecon ).^(2) ), 2 );
    
    % Save the intensity at that radius
    sMap( rCount, : ) = intensity;
    
end

reconTime = toc;
disp([ 'Reonstruction Time: ', num2str(reconTime) ]);

%% Plot pressure at each radius
figure()
hold all;

% Retain only the actual channels of S (i.e., not the pad channels)
startIndex = padLength + 1;
endIndex = length( S0( :, 1 ) ) - padLength;
theta = theta( startIndex : endIndex );
sMapActual = sMap( :, startIndex : endIndex );
sMapNorm = sMapActual./max( abs( sMapActual(:) ) );
sMapDb = 20.*log10( sMapNorm );

rOfInterest = 20E-3; % [m]
rIndex = find( rVector > rOfInterest, 1);
plot( theta.*180./pi, sMapDb( rIndex, : ) );

xlabel( 'Azimuth [deg]' );
ylabel( 'Intensity [AU]' );
ylim( [0.5, 1] );

figure()
hold on;

[rPlot, thetaPlot] = meshgrid( rVector, theta );
zPlot = rPlot.*cos( thetaPlot );
xPlot = rPlot.*sin( thetaPlot );
% pcolor( zPlot.*1E3, xPlot.*1E3, sMapActual' );
pcolor( zPlot.*1E3, xPlot.*1E3, sMapDb' );
shading flat;

%% Get gradient on re-sampled grid

% Resample to make plot less busy
if plotGradient
    factor = 3; % Factor to reduce num points by
    rStep = factor;
    thetaStep = round( ...
        ( length(mVector)./length(rVector)).*factor ...
        );
    
    % If r is descreasing, flip it so we have positive outward
    if rVector(end) < rVector(1)
        rVector = fliplr( rVector );
        % Warn user (this should be the last time we use rVector, but just
        % in case...
        disp( 'WARNING: Flipping rVector!' );
    end    
    
    % Create vectors of indices to keep
    [rDim, thetaDim] = size( sMapActual );
    rKeepIndices = rStep : rStep : rDim - rStep;
    thetaKeepIndices = thetaStep : thetaStep : thetaDim - thetaStep;
    if max( rKeepIndices ) > length( theta)
        thetaKeepIndices(end) = length( theta );
    end
    % Create interpolation points
    [ rInterp, thetaInterp ] = meshgrid( ...
        rVector( rKeepIndices ), theta( thetaKeepIndices ) );
        % Now get plotting vectors
    zInterp = rInterp.*cos(thetaInterp);
    xInterp = rInterp.*sin(thetaInterp);
    % Interpolate function
    sMapResampled = interp2( rPlot, thetaPlot, sMapActual', ...
        rInterp, thetaInterp );
    
    % Get the gradient
    [zComp, xComp] = gradient( sMapResampled );
    
    % Plot!
    quiver( zInterp.*1E3, xInterp.*1E3, zComp, xComp, ...
        'Color', 'w');
end

% Plot source positions
numSources = data.numSources;
sourcePositions = data.sourcePositions;
for sourceCount = 1:numSources
    
   zSource = sourcePositions( sourceCount, 1 ) - z0; 
   xSource = sourcePositions( sourceCount, 2 ) - x0;
   plot( zSource.*1E3, xSource.*1E3, 'ro' );
    
end

axis equal;
xlabel( 'Axial Distance [mm]' );
ylabel( 'Transverse Distance [mm]' );
