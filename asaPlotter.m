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
pchannels = 4*128;
stps = 128;

%% Set parameters
bin = 'n';
deld = 1;
numTargets = data.trgt; % Number of targets
dx = data.dx*1E3; % Element spacing [mm]

% Medium proerties
c = 1482; % [m/s]
rho0 = 998; % [kg/m^3]

% Set location of the probe [mm]
loc = dx*(data.yDim - 5) - data.excit_loc(numTargets,2,1)*dx;

% Compute the time step
dt =data.t(2) - data.t(1); % [s]
t = data.t; % Time vector [s]

% Get the received data for a row of points at the desired x-location
xReceiverPosition = data.excit_loc(numTargets,1,1);
rfdata = squeeze( ...
    data.aedata(xReceiverPosition, :, :) );

% Get the size of the received US data array. This array has dimension:
%  (index of sensor) by (time index)
ss = size(rfdata);
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
for mm = 1:channels
    SNR = 120; % [dB] 120 - Almost no noise, 58/62 - High noise
    rf(mm, :) = addGaussianWhiteNoise(rf1(mm,:),SNR);
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

%% power spectrum estimation
Ndata=length(rf);
fftd = fft( double( rf(channels/2,:)));
pw = fftd.*conj(fftd)/(Ndata.^2);
psd = 2.*pw( 1:floor(Ndata/2) );
Fs = 1/(dt.*1E3);
frqax = (0:Ndata/2-1).*Fs;

subplot(2,1,2)
plot(frqax./1E6, psd)
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

% Set up harminic frequency bins
[~,h1] = min( abs(frqax - 2.8E5) );
[~,h2] = min( abs(frqax - 14.8E5) );
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
p=rfasa';

% Take the FFT of the padded data on each channel
aa1=(fft(p));

% Now for each spatial frequency bin
for mm=1:ss(1)
    
    % Get the center frequency
    fc = frqax(fbins(mm)); % [Hz]
    w = 2*pi*fc; % [rad/s]
    
    % ff the data at fc
    [zz, ff] = min( abs(fk-fc) );
    
    % Get complex data at fc for each channel - CC based on single FFT outside loop
    P = aa1(ff+1,:);
    
    % Convert P into K-space
    P=fftshift(fft(P));
    
    % Image reconstruction algorithm -------------
    dx = x(2)-x(1);
    dk = 2.*pi./dx; % Wavenumber increment
    startValue = ( -floor( length(x)/2 ) );
    endValue = ( ceil( length(x)/2 ) - 1 );
    
    % Create wavenumber vector
    k = (startValue:endValue).*dk./length(x);
    
    % Initialize array to hold final image
    asa = zeros( length(z), pchannels);
    
    %     for lp=1:length(z)
    %         asa(lp,:)=(ifft(ifftshift(P.*exp(1i.*(z(lp)).*sqrt(w.^2./c.^2-k.^2)))))';
    %     end;
    
    % vectorized version
    Pv = repmat(P,length(z),1); %128x512
    kv = repmat(k,length(z),1); %128x512
    zv = repmat(z,length(x),1); %128x512
    
    % Get the wavenumber in the propagation direction
    kz = sqrt( w.^(2)./c.^(2) - kv.^(2) );
    
    % Apply shift to data  and take inverse transform to recover field at
    % source plane
    asa = ifft( ...
        ifftshift( Pv.*exp(1j.*kz.*zv'), 2),[],2 ...
        );
    
    % Get squared magnitude of the signal
    asa = 2*( abs(asa) ).^(2);
    
    % Add contribution of this spatial frequency bin
    asamap = asamap + asa;
end
toc

%% Plot
figure()

asim = imrotate(asamap,-90);
imagesc(-z*1000,x*1000,flip(asim,2));title('AS-PAM','FontWeight','bold')
axis image
% colormap jet;
ylabel('Transeverse [mm]','FontWeight','bold')
xlabel('Axial [mm]','FontWeight','bold')
ylim([-30,30])
caxis([0 max(max(asim))])