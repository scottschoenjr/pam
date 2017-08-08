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
% sourceFile = ...
%     '../data/results/skullData_3bub';
% sourceFile = ...
%     '../data/results/oneBubble/uniform/uniformMedium_1500mps_1bub_4040';
sourceFile = ...
    '../data/results/oneBubble/circularLayer/circularLayer_3000mps_4040_1bub_4045';
% sourceFile = ...
%     '../data/results/oneBubble/stratified/stratified_gaussian_2500pk_005RangeSigma_1bub_4040.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

%% Set parameters

% ------------------------ Other Settings ---------------------------------
fdtdOffset = 19.*data.dx;
padFactor = 1; % 1 - No Pad channels, 2 - Pad by N/2, etc.
constructEntireRegion = 0;
zRange = [50, 70]./1E3 + fdtdOffset; % Set range of reconstruction [m]
isTrulyStratified = 0;
% -------------------------------------------------------------------------

bin = 'y';
deld = 2;
numTargets = data.numSources; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

% Medium proerties
c = 3000; % [m/s]
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
channels = 128;

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
for mm = 1:channels
    SNR = 70; % [dB] 120 - Almost no noise, 58/62 - High noise
    rf(mm, :) = addGaussianWhiteNoise(rf1(mm,:),SNR);
end

% Power spectrum at center channel
centerChannelIndex = round( numSensors./2 );
spectrumAmplitude = abs( fft( arrayData(centerChannelIndex, :) ) );
powerSpectrumNorm = ( spectrumAmplitude./max( spectrumAmplitude ) ).^(2);

Fs = 1./dt;
df = Fs./numTimeSteps;
fVector = 0 : df : ( numTimeSteps - 1 ).*df;
halfwayIndex = floor( numTimeSteps )./2;

%% Compute Field with Angular Spectrum Approach

% Set up number of padding channels (in each direction!)
padChannels = numSensors;
% This is the number of chanels, plus the number of pad channels
totalChannels = round( padFactor.*channels ); 

% Set the number of points at which the field will be calculated
zEvaluationPoints = 256;

% Initialize angular spectrum map
asamap = zeros(zEvaluationPoints, totalChannels);
asamapRef = zeros(zEvaluationPoints, totalChannels);

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

% Get the sound speed vector
dataC = data.c; % As a function of z *from left* (not receiver!)
dataZ = data.xPositionVector; % [m]

% Slice out the portion between the array and the greatest depth
zStartIndex = find( ... % Will fail if z depth too big
    dataZ > (data.recPosition - (z(end) + dx)), 1 ); % Add to make sure zClipped has bigger range than z
zEndIndex = find( ...
    dataZ > (data.recPosition - z(1)), 1 );
% Use profile along any row (since it's stratified). I chose 117 since it's
% more likely to catch eyes than 1.
cClipped = dataC( 117, zStartIndex : zEndIndex );
% %%%%%%% DEBUG %%%%%%%
if ~isTrulyStratified
    cClipped = mean( dataC(:, zStartIndex : zEndIndex ), 1 );
    disp( 'Averaging sound speed in x.' );
end
% %%%%%%%%%%%%%%%%%%%%%%
zClipped = data.recPosition - dataZ( zStartIndex : zEndIndex );

% Interpolate sound speed to evaluated z depths
c_z = interp1( zClipped, cClipped, z );
c = mean( c_z );

% To reconstruct at a subset of points
if ~constructEntireRegion
    
    % Adjust offset
    if ~isTrulyStratified
        zSpanNew = max(zRange) - min(zRange);
        zSpanOld = max(z) - min(z);
        ratio = zSpanNew./zSpanOld;
        fdtdOffset = ratio.*fdtdOffset;
    end
    
    zStartIndex = find( z > zRange(1), 1 );
    zEndIndex = find( z > zRange(2), 1 );
    z = z( zStartIndex : zEndIndex );
    c_z = c_z( zStartIndex : zEndIndex );
    asamap = zeros(length(z), totalChannels);

end

% Compute phase delays for all frequncies and positions
tic;
phiMat = zeros( length(ss(1)), length(z), length(x) );
for fCount = 1 : ss(1)
    
    % Get the center frequency
    fc = fVector(fbins(fCount)); % [Hz]
    omega = 2*pi*fc; % [rad/s]
    
    % Assemble wavenumber vector
    dx = x(2)-x(1);
    dk = 2.*pi./dx; % Wavenumber increment
    startValue = ( -floor( length(x)/2 ) );
    endValue = ( ceil( length(x)/2 ) - 1 );
    
    % Create wavenumber vector
    k = (startValue:endValue).*dk./length(x);
    
    % Compute the propagating wavenumber
    kv = repmat(k, [length(z), 1]);

    % Get the wavenumber in the propagation direction
    kz = sqrt( omega.^(2)./c.^(2) - kv.^(2) );    
    
    % Compute the phase delay
    c0 = c; % Reference sound speed
    cPrime = fliplr( c_z - c0 ); % Measure back from receiver
    mu_z = ( 1 + cPrime./c0 ).^(-2);
    lambda_z = (omega.^(2)./c.^(2)).*( 1 - mu_z );
          
    % Now compute the phase delay for each receiver for each depth
    for zCount = 1:length(z)
        
        integrand = lambda_z(1:zCount).*dz;
        phiMat(fCount, zCount, :) = (1./(2.*kz(zCount, :))).*sum( integrand );
       
    end
    
end
% Set pure imaginary elements to 0
phiMat = mod( real(phiMat), 2.*pi );
clc;
computationTime = toc;
displayString = ...
    sprintf( 'Phase Delay Computation Time: %6.2f s', computationTime );
disp( displayString );

% Time reconstruction
tic

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
    
    % Vectorize so that we can perform the back propagation at each depth z
    % in a single step.
    Pv = repmat(P, [length(z), 1]); % totalChannels x zEvaluationPoints
    kv = repmat(k, [length(z), 1]); % 
    zv = repmat(z, [length(x), 1]); % 

    % Get the wavenumber in the propagation direction
    kz = sqrt( omega.^(2)./c.^(2) - kv.^(2) );     
    
    % Get additional phase delay
    phi = squeeze( phiMat( fCount, :, : ) );
    
    % Add in phase correction
    asa = ifft( ...
        ifftshift( Pv.*exp( 1i.*kz.*zv' + 1i.*phi ), 2), [], 2 ...
        );
    
    % Get squared magnitude of the signal
    asa = 2*( abs(asa) ).^(2);
    
    % Add contribution of this frequency
    asamap = asamap + asa;
    
end

% Display computation time
computationTime = toc;
displayString = ...
    sprintf( 'ASA Method Computation Time: %6.4f s', computationTime );
disp( displayString );

%% Plot angular spectrum result

figure()
hold all;

% Account for FDTD offset
z = z - fdtdOffset;

[ xAsaPlot, zAsaPlot ] = meshgrid( x, z );
pcolor( zAsaPlot.*1E3, xAsaPlot.*1E3, asamap./max( asamap(:) ) );
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

xlabel( 'Distance from Receiver [mm]', 'FontSize', 26 );
ylabel( 'Transverse Distance [mm]', 'FontSize', 26 );

cBarHandle = colorbar;

cBarHandle.Label.String = 'Normalized Intensity';
cBarHandle.Label.FontSize = 22;
cBarHandle.Label.Interpreter = 'latex';
cBarHandle.TickLabelInterpreter = 'latex';

% Get profiles at peak value
[maxValue, maxIndex] = max( asamap(:) );
[maxRow, maxCol] = ind2sub( size(asamap), maxIndex );

% Axial profile
figure()
axialProfile = asamap( :, maxCol );
axialProfileNorm = axialProfile./max(axialProfile);
plot( 1E3.*z, axialProfileNorm, 'k' );
ylim( [0, 1.01] );
xlabel( 'Distance From Receiver [mm]' );
set( gca, 'XDir', 'Reverse' );
ylabel( 'Normalized Intensity' );

fwhm( axialProfileNorm, z.*1E3, 1 );

% Transverse profile
figure()
transProfile = asamap( maxRow, : );
transProfileNorm = transProfile./max(transProfile);
plot( 1E3.*x, transProfileNorm, 'k' );
ylim( [0, 1.01] );
xlabel( 'Transverse Distance [mm]' );
set( gca, 'XDir', 'Reverse' );
ylabel( 'Normalized Intensity' );

% fwhm( transProfileNorm, z.*1E3, 1 );

% Plot phase corrections
run phaseCorrectionPlotter.m