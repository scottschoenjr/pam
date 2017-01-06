%**************************************************************************
%
% Reconstruct Field (2D)
%
%   This script uses a passive ultrasonography algorithm based on back-
%   propagation to reconstruct cavitation from selective images and
%   superimposes them to construct 2D maps of cavitation activity. The
%   spectral response of the cavitation events (temporal information) and
%   of the cavitation map (treatment information) is also extracted.
%
%
% 201008XX - Costas Arvanitis 
% 201607XX - Calum Crake (Frequency Domain Method)
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
    '../data/results/oneBubble/stratifiedMedium_200mps_025mm_1bub.mat';
data = load(sourceFile);
disp(['               ...done (', num2str(toc), ' s).' ] )

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

figure (10)
subplot(1,2,1)
imagesc(1e6*t,0:xtr/(channels-1):xtr,rf)
ylabel('Distance','FontWeight','bold')
xlabel('Time (usec)','FontWeight','bold')
hb = colorbar;
ylabel(hb,'Pressure (Pa)','FontSize',12,'FontWeight','bold')

%% power spectrum estimation
Ndata=length(rf);
fftd=fft(double(rf(channels/2,:)));
pw=fftd.*conj(fftd)/Ndata^2;
psd=2*pw(1:floor(Ndata/2));
df=1/data.t(2);
frqax=(0:Ndata/2-1)*df/Ndata;

subplot(1,2,2)
plot(frqax,psd)
xlim([2e+5,15e+5])
ylabel('Intensity (a.u.)','FontWeight','bold')
xlabel('Normalized Frequency','FontWeight','bold')
title('Power Spectrum','FontWeight','bold')

%% Time domain algorithm
disp('Computing Time Domain Solution');

% z axis (axial)
t_zlength = 1E-3*far_coor;            % image depth (m)
t_z = linspace(0, t_zlength, 128);      % z axis (m) [1x128]

% x axis (transverse)
t_xlength=xtr*0.001;                % x axis length (m) based on size of array
t_x=linspace(-t_xlength/2,t_xlength/2,128);     % x axis (m) [1x128]
vox = xz2vox(t_x*1000,t_z*1000);% for zooming at specific point in the reconstracted image

%vox = xz2vox((lft_el:step:rght_el),near_coor:step:far_coor);% for zooming at specific point in the reconstracted image
del = xtr/channels; % for 64 channles
xs = [ ...
    ( (-del*channels/2):del:0 ), (del:del:(del*(channels/2-1)) ) ...
    ];

% Set ray parameters
ray.fs = (1./dt).*1E-6; % sampling rate (MHz)
ray.x = xs';
ray.N = length(ray.x);
ray.z = zeros(ray.N,1);
ray.y = zeros(ray.N,1);

ray.r = [ray.x, ray.y, ray.z];
ray.A = norm(xs(end)-xs(1));
ray.c=c;

% variables that will be saved
img_frm = zeros(vox.Nz,vox.Nx);
% DC_bias = zeros(vox.Nz,vox.Nx);

fprintf('   ');
zmax = max(vox.z);
rf_ = rf';
[NRF, Nchan] = size(rf'); % in fact this is rf'
tic
% Initialize progress bar
for Ivox = 1:vox.N
    
    % Update status bar every 10 steps
    if ( mod(Ivox, 1000) == 0 )
        fprintf('\b\b\b%3.0f',100*Ivox/vox.N);
    end
    dels = round(sqrt(sum((ones(channels,1) * vox.r(Ivox,:) -ray.r).^2,2)) *1000*ray.fs/ray.c)';
    Ddel = max(dels)-min(dels); % maximum difference in delays
    dels = dels-min(dels) + (0:(Nchan-1))*NRF;
    Nsample=1;
    Ip = 1:Nsample:(NRF-Ddel);
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
    img_frm(Ivox) = sum((sqdels*delrf).^2 - (sqdels*delrf.^2));
end
disp(''); % New line
toc;
%align image to simulation parameters
pus=imrotate(img_frm,-90);
uszz=vox.z;
usxx=vox.x-min(vox.x);
const= data.zDim-(channels*floor(ss(1)/channels));

%% Compute Field with Angular Spectrum Approach
disp('ASA:')

% Set up number of channels
pchannels = 4*128;
stps = 128;

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
    asamap = asamap+asa;
end
toc

%% Frequency domain
% follows same layout as ASA code above, but using algorithm as per
% Haworth et al 2012 'passive imaging with pulsed ultrasound insonations'
%
% variables are prefixed with 'f_' to keep separate from existing code
% new terms follow Haworth 2012:
% f_s = time domain RF data
% f_S = frequency domain RF data
% f_A = apodization term
% f_B = initial image (Haworth 2012 equation 1)
% f_I = image after subtracting DC term (equation 2)
%
% The final map after summing over channels and frequencies is in f_map
%
% calculation over voxels is vectorized by defining grid/element locations
% using complex notation, where j component is z axis. this allows easy
% distance calculation over the whole grid by taking the magnitude of the
% vectors, e.g. abs(f_gridPos-f_elPos)
%
% Calum Crake July 2016


disp('FD:')
selbw = 1; % 1=use bands as per ASA, 2=full BW (slow)

f_channels=128; % use 128 channels (not zero padded)
f_stps=128;     % number of steps in depth
f_rf=rf;        % use original RF (not transposed or zero padded)[128x1787]

% z axis (axial)
f_zlength=1e-3*far_coor;            % image depth (m)
f_z=linspace(0,f_zlength,f_stps);   % z axis (m) [1x128]

% x axis (transverse)
f_xlength=xtr*0.001;             % x axis length (m) based on size of array
f_x=linspace(-f_xlength/2,f_xlength/2,f_channels);     % x axis (m) [1x128]

% combine axes into grid using complex notation, where j component= z axis
f_gridPos = repmat(f_x,[length(f_z) 1]) + 1j*repmat(f_z',[1 length(f_x)]);

switch selbw
    case 1
        disp('Using frequency bins as per ASA')
        % harmonic frequency bins - same as per ASA
        % where frqax is the frequency axis (Hz) from RF data calculated above
        [~,f_h1]=min(abs(frqax-2.8e+005));  % finds closest fft indices
        [~,f_h2]=min(abs(frqax-14.8e+005)); % within frequency axis freqax
        f_fbins=(f_h1:f_h2)';    % list of frequency indices to run [50x1]
        f_ss=size(f_fbins);      % f_ss(1) = number of frequency indices to run (50)
    case 2
        disp('Using full BW as per TD')
        % harmonic frequency bins - full bandwidth
        % this is much slower and doesn't seem to affect the image
        % Haworth paper uses discrete frequency bands so probably the above
        % approach is more appropriate
        f_h1=1;
        f_h2=length(frqax);
        f_fbins=(f_h1:f_h2)';    % list of frequency indices to run [893x1]
        f_ss=size(f_fbins);      % f_ss(1) = number freq indices to run (893)
end

f_B = zeros(size(f_gridPos,1),size(f_gridPos,2),f_ss(1));
f_I = f_B;

f_A = 1; % apodization term - simplest case uses unity term everywhere (see Haworth 2012 eqn 3a)

f_dt = data.t(2); % time step (s)
f_fk =(1/(f_dt))*linspace(0,1,length(f_rf)); % frequency vector [1x1787]

tic
f_s=f_rf'; % rf data for 128 channels- outside loop as per revised ASA code
f_aa1=(fft(f_s)); % fft of whole data  - outside loop to be consistent with revised ASA code
fprintf('   ');
for f_mm=1:f_ss(1)           % loop over frequency bins
    fprintf('\b\b\b%3.0f',100*f_mm/f_ss(1));
    f_fc=frqax(f_fbins(f_mm)); % current frequency (Hz) using bin number in freqax calculated above
    f_w=2*pi*f_fc;           % current angular frequency (rad/s)
    % find the position of fc
    [~, f_ff]=(min(abs(f_fk-f_fc))); % index number where fk is closest to fc
    % eg for f_mm=1, fc=2.9184e+05, f_ff=13 > f_fk(13)=2.9200e+05
    % extract the complex data P from each channel at fc
    f_S=zeros(1,f_channels); % initialize variable to hold current freq domain data [1x128]
    f_m1 = zeros(size(f_gridPos)); % used to sum map for current freq
    f_m2 = 0;                      % second term Haworth 2012 equation 2
    % loop over channels - Haworth 2012 equation 1
    for f_gg=1:f_channels
        % extract the data at fc for current element (one complex number)
        f_S(f_gg)=f_aa1(f_ff+1,f_gg);
        % get current element loc in complex notation (j component = z axis)
        f_elPos = ray.x(f_gg,1)/1000 + 1j*ray.z(f_gg,1)/1000; % element pos (m)
        
        % Image formation algorithm: f_m1 is the summation term in Haworth
        % 2012 equation 1. The abs()^2 operation is taken below after
        % summing over the elements in this loop
        f_m1 = f_m1 + (f_S(f_gg).*f_A.*exp(1j*f_w*abs(f_gridPos-f_elPos)/c));
        
        % f_m2 is the summation term in Haworth 2012 equation 2. It is
        % subtracted from the image produced by equation 1 below after
        % summing over the elements in this loop
        f_m2 = f_m2 + abs(f_S(f_gg).*f_A).^2; % second sum term in Haworth 2012 equation 2
        
        % uncomment below for nice (but slow!) visualization
        % figure(99); imagesc(abs(f_m1)); title(['freq ' num2str(f_mm) ', ch ' num2str(f_gg)]); pause(0.01);
    end
    f_B(:,:,f_mm) = abs(f_m1).^2;         % Haworth 2012 equation 1
    f_I(:,:,f_mm) = f_B(:,:,f_mm) - f_m2; % Haworth 2012 equation 2
end
% finally, sum over frequencies for the final image as per Haworth 2012 (see text below equation 2)
f_map = sum(f_I,3)./f_channels; % normalized for number of channels
toc
%% final images for the different pam algortihms
fim=imrotate(f_map,90);
asim=imrotate(asamap,-90);
tdim=pus;

plotfwhm = 0;
%% stats

% Time domain stats
[~, tdloc1]=max(max(tdim));
[~, tdloc2]=max(tdim(:,tdloc1));
[~, tdloc3]=max(max(tdim'));
[tdval, tdloc4]=max(tdim(tdloc3,:));
[mintdval1, ~]=min(tdim(:,tdloc1));
[mintdval2, ~]=min(tdim(tdloc3,:));
% Spatial resolution analysis: FWHM estimation
tdpam_ax=(tdim(tdloc3,:)-mintdval2)/(tdval-mintdval2);
tdpam_tr=(tdim(:,tdloc1)-mintdval1)/(tdval-mintdval1);
%FWHM
tdwdth=vox.z(tdpam_ax-0.49>=0);
td_fwhm=tdwdth(end)-tdwdth(1);
%SNR
snr_tr_tdim=1/std(tdpam_tr(tdloc2-12*ceil(1/step):tdloc2-4*ceil(1/step)));

td_fwhm_ax = fwhm(tdpam_ax,vox.z,1); % show one plot as example
td_fwhm_tr = fwhm(tdpam_tr,vox.x,plotfwhm);

% AS stats
[~, asloc1]=max(max(asim));
[~, asloc2]=max(asim(:,asloc1));
[~, asloc3]=max(max(asim'));
[asval, ashloc4]=max(asim(asloc3,:));
[minasval1, ~]=min(asim(pchannels/2-channels/2:pchannels/2+channels/2,asloc1));
[minasval2, ~]=min(asim(asloc3,:));
%td_fwhm=0
aspam_ax=(asim(asloc3,:)-minasval2)/(asval-minasval2);
aspam_tr=(asim(:,asloc1)-minasval1)/(asval-minasval1);
aswdth=z(aspam_ax-0.49>=0);
as_fwhm=(aswdth(end)-aswdth(1))*1000;
% Image Quality metrics: SNR
trnsvprof_asim=asim(:,ashloc4)/asval;
snr_tr_asim=1/std(aspam_tr(asloc2-12*ceil(1/step):asloc2-4*ceil(1/step)));

as_fwhm_ax = fwhm(aspam_ax,z*1000,plotfwhm);
as_fwhm_tr = fwhm(aspam_tr,x*1000,plotfwhm);

% FD stats
[~, fdloc1]=max(max(fim));
[~, fdloc2]=max(fim(:,fdloc1));
[~, fdloc3]=max(max(fim'));
[fdval, fdloc4]=max(fim(fdloc3,:));
[minfdval1, ~]=min(fim(:,fdloc1));
[minfdval2, ~]=min(fim(fdloc3,:));
% Spatial resolution analysis: FWHM estimation
fdpam_ax=(fim(fdloc3,:)-minfdval2)/(fdval-minfdval2);
fdpam_tr=(fim(:,fdloc1)-minfdval1)/(fdval-minfdval1);
%FWHM
fdwdth=vox.z(fdpam_ax-0.49>=0);
fd_fwhm=fdwdth(end)-fdwdth(1);
%SNR
snr_tr_fim=1/std(fdpam_tr(fdloc2-12*ceil(1/step):fdloc2-4*ceil(1/step)));

fd_fwhm_ax = fwhm(fdpam_ax,vox.z,plotfwhm);
fd_fwhm_tr = fwhm(fdpam_tr,vox.x,plotfwhm);

table_fwhm=[td_fwhm;fd_fwhm;as_fwhm]

table_fwhm_new_ax=[td_fwhm_ax;fd_fwhm_ax;as_fwhm_ax]
table_fwhm_new_tr=[td_fwhm_tr;fd_fwhm_tr;as_fwhm_tr]

table_snr=[snr_tr_tdim;snr_tr_fim;snr_tr_asim]

%% Disp images and data
fig=figure(12);
set(fig, 'position', [25 50 1900 900])

% figure (12)
subplot(2,4,1)
imagesc(1e6*t,1:128,rf)
axis image
ylabel('Channel','FontWeight','bold')
xlabel('Time [\mus]','FontWeight','bold')
title('RF data','FontWeight','bold')

subplot(2,4,5)
plot(frqax./1E6,psd./max(psd),'k')
xlim([0.2,1.5])
ylabel('Normalized Intensity','FontWeight','bold')
xlabel('Frequency [MHz]','FontWeight','bold')
title('Power Spectrum','FontWeight','bold')

subplot(2,4,2)
imagesc(loc-z*1000,x*1000,flip(asim,2));title('AS-PAM','FontWeight','bold')
axis image
% colormap jet;
ylabel('Transeverse [mm]','FontWeight','bold')
xlabel('Axial [mm]','FontWeight','bold')
ylim([-30,30])
caxis([0 max(max(asim))])

subplot(2,4,3) % CC added fourth column for FD-PAM plot
imagesc(loc-f_z*1000,f_x*1000,flip(fim,1));
title('FD-PAM','FontWeight','bold')
axis image
% colormap jet;
ylabel('Transverse [mm]','FontWeight','bold')
xlabel('Axial [mm]','FontWeight','bold')
ylim([-30,30])
caxis([0 max(max(fim))])

subplot(2,4,4)
imagesc((loc-uszz),usxx-(data.excit_loc(numTargets,3,1))*dx,flip(tdim,2))
axis image
% colormap gray;
ylabel('Transeverse [mm]','FontWeight','bold')
xlabel('Axial [mm]','FontWeight','bold')
title('TD-PAM','FontWeight','bold')
ylim([-30,30])
caxis([0 max(max(tdim))])

subplot(2,4,6)
plot (loc-z*1000-(z(2)-z(1))*1000,flip(aspam_ax,2),'r')
hold on
plot (loc-f_z*1000-(f_z(2)-f_z(1))*1000,fdpam_ax,'k:'); % CC added FD line
plot ((loc-uszz)+uszz(3),tdpam_ax,'k');
title('Axial Line Profile','FontWeight','bold')
xlabel('Distance From Source [mm]','FontWeight','bold')
ylabel('Normalized Intensity','FontWeight','bold')
xlim([-30,30])
ylim([0,1])
% legend(['AS-PAM';'FD-PAM';'TD-PAM']) % CC added FD line
hold off

subplot(2,4,7:8)
plot (x*1000-(x(2)-x(1))*1000,aspam_tr,'r')
hold on
plot (f_x*1000,flip(fdpam_tr),'k:'); % CC added FD line
plot (usxx-(data.excit_loc(numTargets,3,1))*dx+usxx(2),tdpam_tr,'k');
title('Transverse Line Profile','FontWeight','bold')
legend(' Angular Spectrum',' Frequency Domain',' Time Domain', ...
    'Location', 'SouthEastOutside' ) % CC added FD line
xlabel('Distance From Source [mm]','FontWeight','bold')
ylabel('Normalized Intensity','FontWeight','bold')
xlim([-30,30])
ylim([0,1])
hold off
