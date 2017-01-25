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
    '../data/results/oneBubble/layeredMedium_3000mps_04mm_1bub_30mmoffset';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

%% Set parameters

% Set sound speed method to use
% 0 - Nothing, just use uniform c
% 1 - Account for layer (see that code to adjust parameters
% 2 - Use Averaged (in z) sound speed
% 3 - Use stratified medium result (once it works...)
soundSpeedMethod = 2;

bin = 'y';
deld = 2;
numTargets = data.numSources; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

% Medium proerties
c = 1500; % [m/s]
rho0 = 1000; % [kg/m^3]

% Set location of the probe [mm]
loc = dx*(data.yDim - 5) - data.excit_loc(numTargets,2,1)*dx;

% Compute the time step
t = data.t; % Time vector [s]
dt = data.t(2) - data.t(1); % Time step [s]

% Get the x-positions of the sensors
y = data.yPositionVector; % [m]

% Get the received data
arrayData = data.aedata;

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

%% Add noise to the rf data and display them
rf = rf1;
% for mm = 1:channels
%     SNR = 120; % [dB] 120 - Almost no noise, 58/62 - High noise
%     rf(mm, :) = addGaussianWhiteNoise(rf1(mm,:),SNR);
% end

% Plot raw data and power spectrum
figure (10)

subplot(1,2,1)
[tColorPlot, yColorPlot] = meshgrid( t, y );
pcolor( tColorPlot.*1E6, yColorPlot.*1E3, arrayData );
shading flat;

ylabel('Position [mm]' )
xlabel('Time [$\mu$s]')
cdHandle = colorbar;
ylabel(cdHandle, 'Pressure [Pa]', 'Interpreter', 'latex')

% Power spectrum at center channel
centerChannelIndex = round( numSensors./2 );
spectrumAmplitude = abs( fft( arrayData(centerChannelIndex, :) ) );
powerSpectrumNorm = ( spectrumAmplitude./max( spectrumAmplitude ) ).^(2);

Fs = 1./dt;
df = Fs./numTimeSteps;
fVector = 0 : df : ( numTimeSteps - 1 ).*df;
halfwayIndex = floor( numTimeSteps )./2;

subplot(1,2,2)
plot(fVector(1:halfwayIndex)./1E6, powerSpectrumNorm(1:halfwayIndex), 'k' );
ylim([0, 1]);
xlim([0, 5]);
ylabel('Normalized Power','FontWeight','bold');
xlabel('Frequency [MHz]','FontWeight','bold');

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

% Now for each frequency bin
for fCount = 1:ss(1)
    
    % Get the center frequency
    fc = fVector(fbins(fCount)); % [Hz]
    omega = 2*pi*fc; % [rad/s]
    
    % If the averaged sound speed is to be used, calculate it here
    if soundSpeedMethod == 2
        
        % Get simulation data
        dataC = data.c; % Total field
        sourcePositions = data.sourcePositions; % Position of sources
        zStart = data.recPosition - z(end); % Position of leftmost
        zEnd = data.recPosition - z(1); % Position of receiver
        
        % Get z-indidces of data to keep
        dataZ = data.xPositionVector;
        zStartIndex = find( dataZ > zStart, 1 );
        zEndIndex = find( dataZ > zEnd, 1 );
        
        % Get portion of cField and average it
        cField = dataC( :, zStartIndex : zEndIndex );
        % Average sound speed in each plane
        cFieldMeanPlane = mean( cField, 1 ); 
        % Average sound speed in each plane, interpolated to the
        % appropriate length
        numPointsRaw = zEndIndex - zStartIndex + 1; % Points in data
        rawPoints = linspace(0, 1, numPointsRaw);
        numPointsInterp = length(z); % Number of points to interpolate to
        interpPoints = linspace(0, 1, numPointsInterp );
        cFieldMeanPlaneInterp = interp1( ...
            rawPoints, cFieldMeanPlane, interpPoints );  
        
        % New averaging...
        % cProfile contains the average sound speed from the first to the
        % current entry [e.g., cProfile(10) = mean(cFieldMeanPlane(1:10))]
        cProfile = meanOfPreviousEntries( cFieldMeanPlaneInterp );
        % But since we're measuring from back from the receiver array, we
        % have to flip the profile lengthwise
        cProfile = fliplr( cProfile );
        % Repeat the profile along x for ASA reconstruction
        c = ( repmat(cProfile, [length(x), 1]) )';
        
        cMean = mean( cFieldMeanPlane ); % Get the average of all planes

        % Set sound speed to this average
        c = cMean;
        
    end
    
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
    
    % Initialize array to hold final image
    asa = zeros( length(z), totalChannels);
    
    %     for lp=1:length(z)
    %         asa(lp,:)=(ifft(ifftshift(P.*exp(1i.*(z(lp)).*sqrt(w.^2./c.^2-k.^2)))))';
    %     end;
    
    % Vectorize so that we can perform the back propagation at each depth z
    % in a single step.
    Pv = repmat(P, [length(z), 1]); % totalChannels x zEvaluationPoints
    kv = repmat(k, [length(z), 1]); % 
    zv = repmat(z, [length(x), 1]); % 

    % Get the wavenumber in the propagation direction
    kz = sqrt( omega.^(2)./c.^(2) - kv.^(2) );
           
    % Apply shift to data  and take inverse transform to recover field at
    % each z evaluation point.
    asa = ifft( ...
        ifftshift( Pv.*exp(1j.*kz.*zv'), 2), [], 2 ...
        );
    
    % If we're propagating to single layer...
    if soundSpeedMethod == 1
       
        % Layer Properties ---------------------
        layerStart = 32; % [mm], from receiver;
        layerThickness = 4; % [mm]
        layerC = 3000; % [m/s]
        %----------------------------------------
        
        % Get indices for each region
        zmm = z.*1E3;
        inds1 = find( zmm < layerStart );
        inds2 = find( ...
            zmm >= layerStart & ...
            zmm < (layerStart + layerThickness) ...
            );
        inds3 = find( zmm >= (layerStart + layerThickness) );
        
        % Get position vectors
        z1 = z( inds1 );
        zv1 = repmat(z1, [length(x), 1]); % 
        z2 = z( inds2 ) - z( max(inds1) );
        zv2 = repmat(z2, [length(x), 1]); % 
        z3 = z( inds3 ) - z( max(inds2) );
        zv3 = repmat(z3, [length(x), 1]); % 
        
        % Get wavenumber vectors
        kv1 = repmat(k, [length(z1), 1]);
        kz1 = sqrt( omega.^(2)./c.^(2) - kv1.^(2) );
        
        kv2 = repmat(k, [length(z2), 1]);
        kz2 = sqrt( omega.^(2)./layerC.^(2) - kv2.^(2) ); % Use layer c
        
        kv3 = repmat(k, [length(z3), 1]);
        kz3 = sqrt( omega.^(2)./c.^(2) - kv3.^(2) );
        
        % Now get angular spectrum of pressure at receiver array
        AS0 = P;
        AS0v = repmat(AS0, [length(z1), 1]);
        
        % Porpagate to layer
        AS1 = AS0v.*exp( 1j.*kz1.*zv1' );
        
        % Get field at layer and initialize the pressure field in the layer
        AS2v = repmat(AS1(end, :), [length(z2), 1]);
        
        % Propagate to other side of layer
        AS2 = AS2v.*exp( 1j.*kz2.*zv2' );
        
        % Repeat on other side of layer
        AS3v = repmat(AS2(end, :), [length(z3), 1]);
        
        % Propagate to end of domain
        AS3 = AS3v.*exp( 1j.*kz3.*zv3' );
        
        % Now concatenate results       
        fullField = [ AS1 ; AS2 ; AS3 ]; % Concatenate
        asa = ifft( ...
            ifftshift( fullField, 2), [], 2 ...
            );
        
    end
    
    % If we have a stratified medium, compute phase delay to add
    if soundSpeedMethod == 3;
        
        % Get sound speed field and take a slice of it along z (any index 
        % for x is fine since stratified
        cField = data.c;
        c_z = cField( centerChannelIndex, : );
        
        % Interpolate to depth computation points
        simZ = data.xPositionVector;
        c_z = interp1( simZ, c_z, z );
        
        % Get auxiliaty functions for that slice
        c0 = 1500;
        k0 = omega./c0;
        mu_z = ( 1 - (c_z - c0)./c0 ).^(-2);
        lambda_z = k0.^(2).*(1 - mu_z);
        
        % Perform integration to each z position
        phi = zeros( zEvaluationPoints, totalChannels );
        for zCount = 1:zEvaluationPoints
            
            zIndex = length(z) - zCount + 1;
            
            % Integrate from 1 to z
            integrand = lambda_z( zCount : end ).*dz;
            phi( zIndex, : ) = (0.5).*(1./kz(zIndex, :))* ...
                sum( integrand )';
            
%             %%%%%%% DENUG %%%%%%%%
%             figure(9988)
%             cla;
%             hold all;
%             plot( z, lambda_z, 'k' );
%             plot( z(zIndex), lambda_z(zIndex), 'ro' );
%             ylabel( '$$\lambda(z)$$' );
%             xlabel( '$$z$$' );
%             drawnow;
%             %%%%%%%%%%%%%%%%%%%%%%%
        end
        
        %%%%%%% DEBUG %%%%%%%
        figure(9989);
        subplot( 2, 2, 1 )
        pcolor( mod( real(phi), 2.*pi ) );
        caxis( [0, 2.*pi] );
        colorbar;
        shading flat;
        title( sprintf( '%3d kHz', floor(fc./1E3) ) );
         
        subplot( 2, 2, 3 )
        pcolor( real(1./kz) );
        caxis( [0, 1200] );
        colorbar;
        shading flat;
        title( sprintf( '%3d kHz', floor(fc./1E3) ) );
        
        drawnow;
        %%%%%%%%%%%%%%%%%%%%%
        
        % Add in phase shift to result
        asa = ifft( ...
            ifftshift( Pv.*exp( 1j.*kz.*zv' ).*exp( -phi ), 2), [], 2 ...
            );
        
    end
    
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

%%%%%%%% DEBUG %%%%%%%%%%%
if soundSpeedMethod == 4;
    figure(999)
    subplot( 2, 1, 1 )
    plot( z.*1E3, c_z, 'k' );
    ylabel( '$c$ [m/s]' );
    subplot( 2, 1, 2 )
    plot( z.*1E3, mu_z, 'k' );
    ylabel( '$\mu$' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot angular spectrum result

figure()
hold all;

[ xAsaPlot, zAsaPlot ] = meshgrid( x, z );
pcolor( zAsaPlot.*1E3, xAsaPlot.*1E3, asamap );
shading flat;

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

cBarHandle = colorbar;

cBarHandle.Label.String = '';
cBarHandle.Label.FontSize = 16;
cBarHandle.Label.Interpreter = 'latex';
cBarHandle.TickLabelInterpreter = 'latex';

% Get axial profile
figure()
hold all
middleIndex = find( x > 0, 1 );
maxValue = max( max( asamap( :,  middleIndex-5:middleIndex + 5 ) ) );
for cpCount = middleIndex-5:middleIndex + 5
    centerProfile = asamap( :, cpCount )./maxValue;
    plot( z, centerProfile + (cpCount - middleIndex).*1.5, 'k' );
end
set(gca, 'XDir', 'Reverse' )
% centerProfileNorm = centerProfile./max(max(centerProfile));
% plot( z, centerProfileNorm, 'k' );
% ylim( [0, 1.01] );
