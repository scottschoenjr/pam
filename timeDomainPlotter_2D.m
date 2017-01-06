%**************************************************************************
%
% Time Domain Field Reconstruction (2D)
%
%   This script uses a passive acoustic mapping (PAM) based on back-
%   propagation to reconstruct cavitation from selective the data collected
%   by the simulated line array. 
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
    '../data/results/oneBubble/layeredMedium_3000mps_04mm_3bub_rec75mm';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

%% Set parameters
bin = 'n';
deld = 1;
numTargets = data.numSources; % Number of sources
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

% Get the received data for a row of points at the desired x-location
xReceiverPosition = data.excit_loc(numTargets,1,1);
% rfdata = squeeze( ...
%     data.aedata(xReceiverPosition, :, :) );
arrayData = data.aedata;

% Get the size of the received US data array. This array has dimension:
%  (index of sensor) by (time index)
ss = size(arrayData);
numSensors = ss(1);
numTimeSteps = ss(2);

% Set the number of receiving channels
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

figure (10)

subplot(1,2,1)
imagesc(1e6*t,0:xtr/(channels-1):xtr,rf)
ylabel('Distance','FontWeight','bold')
xlabel('Time (usec)','FontWeight','bold')
hb = colorbar;
ylabel(hb,'Pressure (Pa)','FontSize',12,'FontWeight','bold')

%% power spectrum estimation
% Ndata=length(rf);
% fftd=fft(double(rf(channels/2,:)));
% pw=fftd.*conj(fftd)/Ndata^2;
% psd=2*pw(1:floor(Ndata/2));
% df=1/data.t(2);
% frqax=(0:Ndata/2-1)*df/Ndata;
% 
% subplot(1,2,2)
% plot(frqax,psd)
% xlim([2e+5,15e+5])
% ylabel('Intensity (a.u.)','FontWeight','bold')
% xlabel('Normalized Frequency','FontWeight','bold')
% title('Power Spectrum','FontWeight','bold')

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

%% Time domain algorithm
disp('Computing Time Domain Solution');

% Set the number of points at which the field will be calculated
zEvaluationPoints = 128;

% Set the image depth
zlength = 1E-3*far_coor; % [m]

% Create array of image depths
z = linspace(0, zlength, zEvaluationPoints); % [m]
dz = z(2) - z(1);

% x axis (transverse)
xlength = xtr*1E-3;   % x axis length (m) based on size of array
x = linspace(-xlength/2, xlength/2, channels);     % x axis (m) [1x128]

% Get voxel object "for zooming at specific point in the reconstracted
% image"
vox = xz2vox(x*1E3, z*1E3);

del = xtr/channels; % for 64 channles
xs = [ ...
    ( (-del*channels/2):del:0 ), (del:del:(del*(channels/2-1)) ) ...
    ];

% Set ray parameters
ray.fs = (1./dt).*1E-6; % Sampling rate [MHz]
ray.x = xs';
ray.N = length(ray.x);
ray.z = zeros(ray.N,1);
ray.y = zeros(ray.N,1);

ray.r = [ray.x, ray.y, ray.z];
ray.A = norm(xs(end)-xs(1));
ray.c = c;

% Initialize time domain reconstructed image
img_frm = zeros(vox.Nz, vox.Nx);

fprintf('   ');
zmax = max(vox.z);
rf_ = rf';
[NRF, Nchan] = size(rf'); % in fact this is rf'
numPixels = vox.N;

% Time computation
tic

% Compute intensity for each pixel
for pixelCount = 1:numPixels
    
    % Update status bar every 200 steps
    if ( mod(pixelCount, 200) == 0 )
        clc;
        fractionDone = pixelCount./numPixels;
        timeTaken = toc;
        timeLeft = (1./fractionDone - 1).*timeTaken;
        fprintf('%3.0f Percent Done, about %4.2f s Left', ...
            100*fractionDone, timeLeft );
    end
    
    % Get vectors from each array element to the pixel of interest
    rn = 1E3.*ones(channels, 1)*vox.r(pixelCount, :);
    r = 1E3.*ray.r;
    
    % Get the delay associated with each element for that pixel
    dels = round( sqrt( ...
        sum( (rn - r).^2, 2 ) ...
        )*ray.fs/ray.c ... % Why *ray.fs?
        )';
    
    % Get the maximum difference in delays
    delayRange = max(dels) - min(dels);
    
    dels = dels-min(dels) + (0:(Nchan-1))*NRF;
    Nsample=1;
    Ip = 1:Nsample:(NRF-delayRange);
    NIp = length(Ip);
    delrf = rf_( ones(Nchan,1)*Ip + dels' * ones(1,NIp) );
    
    % If we assume circular spreading, we can add gain based on distance
    addGain = 0;
    if addGain
        sqdels = sqrt(dels);
    else
        % No gain
        sqdels = ones(length(dels),1)';
    end
    
    % Get value for that pixel based on delays
    img_frm(pixelCount) = sum( (sqdels*delrf).^2 - (sqdels*delrf.^2) );
    
end

% Display computation time
clc;
computationTime = toc;
displayString = ...
    sprintf( 'Time Domain Computation Time: %6.2f s', computationTime );
disp( displayString );

%% Align image to simulation parameters and plot results
figure()
hold all;

tdReconstruction = img_frm;
xVector = vox.x - min(vox.x);
zVector = vox.z;
[x, z] = meshgrid( zVector, xVector );

% Plot source positions
zSourcesPlot = data.recPosition - data.sourcePositions( :, 1 );
xSourcesPlot = data.sourcePositions( :, 2 ); % + 10.*dx./1E3;

% Since ASA plotter treats array center as 0 for x-position, shift all our
% x-values by half of the x-span
xSpan = max(data.yPositionVector) - min(data.yPositionVector); % [mm]

% Now issue plotting commands
pcolor( z, x - 1E3.*xSpan./2, tdReconstruction );
shading flat;
plot( zSourcesPlot.*1E3, (xSourcesPlot - xSpan./2).*1E3, 'ro' );

% Format Plot
set( gca, 'XDir', 'Reverse' );
ylim( [-40, 40] );
xlim( [0, 80] );

xlabel( 'Distance from Receiver [mm]' );
ylabel( 'Transverse Distance [mm]' );

cBarHandle = colorbar;

cBarHandle.Label.String = '';
cBarHandle.Label.FontSize = 16;
cBarHandle.Label.Interpreter = 'latex';
cBarHandle.TickLabelInterpreter = 'latex';


