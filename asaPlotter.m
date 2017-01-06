%**************************************************************************
%
% ASA Plotter
%
%   Script loads data from the FDTD simulation and plots the reconstruction
%   of the field using the angular spectrum approach.
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
    '../data/waterdata2lambda_5Bubbles_z.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

% Set up number of channels
pchannels = 16*128; % Number of padding channels
stps = 128;

%% Set parameters
bin = 'y';
deld = 1;
numTargets = data.trgt; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

% Medium proerties
c = 1482; % [m/s]
rho0 = 998; % [kg/m^3]

% Set location of the probe [mm]
offset = data.excit_loc(numTargets,2,1)*dx; 
loc = dx*(data.yDim - 5) - offset;

% Compute the time step
dt = data.t(2) - data.t(1); % [s]
t = data.t; % Time vector [s]

% Get the received data for a row of points at the desired x-location
xReceiverPosition = data.excit_loc(numTargets,1,1);
xReceierPosition = 0;
rfdata = squeeze( ...
    data.aedata(xReceiverPosition, :, :) );

% Get the size of the received US data array. This array has dimension:
%  (index of sensor) by (time index)
ss = size(rfdata);
numSensors = ss(1);
numTimeSteps = ss(2);

% Set the number of receiving channels
channels =160;

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
rf128 = rfdata(:, rcn_x1:end );

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
rf = 0.0*rf1;
for fCount = 1:channels
    SNR = 120; % [dB] 120 - Almost no noise, 58/62 - High noise
    rf(fCount, :) = addGaussianWhiteNoise(rf1(fCount,:),SNR);
end

% Plot received data
figure(10)
subplot(2, 1, 1)

% Get plotting vectors
xStep = xtr/(channels-1);
xMax = xtr;
distance = 0:xStep:xMax;
tStartIndex = rcn_x1;
tEndIndex = rcn_x1 + size( rf, 2 ) - 1;
time_ms = 1E6*t(tStartIndex:tEndIndex);
nSkip = 8;
channelsToPlot = 1:nSkip:size(rf, 1);

% Plot
waterfallPlotHandle = ...
    waterfall(time_ms, distance(channelsToPlot),  rf(channelsToPlot, :));
ylabel('Distance [mm]', 'FontSize', 14);
xlabel('Time [\mus]','FontSize', 14);
% zlabel( 'Pressure [Pa]','FontSize', 14);
colormap([0, 0, 0]);
set( gca, ...
    'ZTick', [], ...
    'CameraPosition', [18, -300, 0.045] ...
    );

%% power spectrum estimation
Ndata=length(rf);
fftd = fft( double( rf(channels/2,:)));
pw = fftd.*conj(fftd)/(Ndata.^2);
psd = 2.*pw( 1:floor(Ndata/2) );
Fs = 1/(dt.*1E3);
frqax = (0:Ndata/2-1).*Fs;

subplot(2,1,2)
plot(frqax./1E6, psd, 'k')
xlim([0,2])
ylabel('Intensity [AU]', 'FontSize', 14)
xlabel('Frequency [MHz]', 'FontSize', 14)
title('Power Spectrum', 'FontWeight','bold')


%% Reconstruction Algorithm

% Initialize map
asamap = zeros(stps, pchannels);
rfasa = padarray( double(rf), pchannels/2 - channels/2 );

% Set the image depth
zlength = 1E-3*far_coor; % [m]
% Create array of image depths (does not have to be same as num. channels)
z = linspace(0, zlength, stps);

% xlength1=data.zDim*data.dx;%array length
xlength1 = xtr*1E-3;
xlength = pchannels*xlength1/channels;
x = linspace( -xlength/2, xlength/2, pchannels );%for frequency sweep video

% Set up harmonic frequency bins
[~, minIndex] = min( abs(frqax - 2.8E5) );
[~, maxIndex] = min( abs(frqax - 14.8E5) );
fbins = (minIndex:maxIndex)';
ss = size(fbins);

% Get time step and vector
dt = data.t(2) - data.t(1); % Time step [s]
t = data.t; % Time vector [s]

% Create frequency vector
Fs = 1./dt; % Sampling frequency [Hz]
fk = Fs.*linspace(0, 1, length(rf));

% Time ASA computation
tic

% Get the padded RF data into the appropriate dimension
p=rfasa';

% Take the FFT of the padded data on each channel
pTilde=(fft(p));

% Now for each frequency bin
for fCount = 1:ss(1)
    
    % Get the center frequency
    fc = frqax( fbins(fCount) ); % [Hz]
    omega = 2*pi*fc; % [rad/s]
    
    % ff the data at fc
    [zz, ff] = min( abs(fk-fc) );
    
    % Get complex data at fc for each channel - CC based on single FFT outside loop
    pTilde_loop = pTilde(ff+1, :);
    
    % Compute angular spectrum
    A_pTilde = fftshift(fft(pTilde_loop));
    
    % Image reconstruction algorithm -------------
    dx = x(2) - x(1);
    dk = 2.*pi./dx; % Wavenumber increment
    startValue = ( -floor( length(x)/2 ) );
    endValue = ( ceil( length(x)/2 ) - 1 );
    
    % Create wavenumber vector
    k = (startValue:endValue).*dk./length(x);
    
    % Initialize array to hold contribution to image for this bin
    asa = zeros( length(z), pchannels);
    
    %     for lp=1:length(z)
    %         asa(lp,:)=(ifft(ifftshift(P.*exp(1i.*(z(lp)).*sqrt(w.^2./c.^2-k.^2)))))';
    %     end;
    
    % vectorized version
    Pv = repmat(A_pTilde, length(z), 1); %128x512
    kv = repmat(k, length(z), 1); %128x512
    zv = repmat(z, length(x), 1); %128x512
    
    %%%%%%%% DEBUG %%%%%%%%
    alpha = 0.1; % [Np/m]
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Get the wavenumber in the propagation direction
    kz = sqrt( (omega.^(2)./c^(2)) - kv.^(2) ) - 1j.*alpha;
    
    % Apply shift to data  and take inverse transform to recover field at
    % source plane
    asa = ifft( ...
        ifftshift( Pv.*exp(1j.*kz.*zv'), 2),[],2 ...
        );
    
    % Get squared magnitude of the signal
    asa = 2*( abs(asa) ).^(2);
    
    % Add contribution of this spatial frequency bin
    asamap = asamap + asa;
    
    %%%%%%%% DEBUG %%%%%%%
%     figure(998)
%     set(gcf, 'Position', [50, 100, 500, 400] );
%     pch = pcolor( real( sin( );
%     set( pch, 'EdgeColor', 'none' );
%     colorbar;
    
    figure(999)
    set(gcf, 'Position', [600, 100, 500, 400] );
    pch = pcolor(asamap);
    set( pch, 'EdgeColor', 'none' );
    pause(0.15);
    %%%%%%%%%%%%%%%%%%%%%
    
end
toc

%% Plot

% Get angular spectrum map to plot
asim = imrotate(asamap,-90);
asaPlotData = flip(asim,2);

axialSlicePos = 34; %34; % [mm]
transSlicePos = -0.5; % [mm]

axialLimits = [-60, 0]; % [mm]
transLimits = [-30, 30]; % [mm]

% Get slices to plot cross-sections
axialIndex = find( z*1000 > axialSlicePos, 1);
transIndex = find( x*1000 > transSlicePos, 1);

axialSlice = asaPlotData( transIndex, :  );
% axialSlice = sum( asaPlotData, 1 );
axialSliceNorm = axialSlice./max(abs(axialSlice));
transSlice = asaPlotData( :, axialIndex );
% transSlice = sum( asaPlotData, 2 );
transSliceNorm = transSlice./max(abs(transSlice));

% Get vectors to plot positions of cross-section
xAxialSlice = [ -z(axialIndex), -z(axialIndex) ]*1000;
yAxialSlice = [ -1E6, 1E6 ];
xTransSlice = [ -1E6, 1E6 ];
yTransSlice = [ x(transIndex), x(transIndex) ]*1000;

% Create plot
figure()
set(gcf, 'Position', [50, 50, 1100, 850] );
hold all;

% Create axis for 2D plot
ax2D = gca;
set( ax2D, 'Position', [0.02, 0.3, 0.7, 0.65] );

imagesc(-z*1000, x*1000, asaPlotData);
plot(xAxialSlice, yAxialSlice, '--w');
plot(xTransSlice, yTransSlice, '--w');
set( ax2D, 'XTickLabel', '' );

axis image
title('Angular Spectrum Reconstruction', 'FontSize', 18)
ylabel('Transeverse [mm]', 'FontSize', 16);
xlim(axialLimits);
ylim(transLimits);
caxis([0 max(max(asim))]);
% 
% % Plot bubble positions
% for sourceCount = 1:numTargets
%     xIndex = data.excit_loc( sourceCount, 2 );
%     xPos = x( xIndex );
%     zIndex = data.excit_loc( sourceCount, 3 );
%     zPos = -z( zIndex );
%     plot(zPos, xPos, 'ro' );
% end

% To the right of colorplot
transverseAxes = axes( 'Position', [0.65, 0.3, 0.21, 0.65] );
plot(transSliceNorm, x*1000, 'k');
xlim([0, 1]);
ylim(transLimits);
set(gca, 'YAxisLocation', 'right' );

% Below colorplot
axialAxes = axes( 'Position', [0.12, 0.08, 0.5, 0.2] );
plot(-z*1000, axialSliceNorm, 'k');
xlim(axialLimits);
ylim([0, 1]);
xlabel('Axial Position [mm]', 'FontSize', 16);