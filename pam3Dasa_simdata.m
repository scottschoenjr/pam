%% from Costas Arvanitis (C) August 2013 BWH
%3D passive ultra-sonography algorithm based on Angular Spectrum Approach
close all;
clear all;
%% parameters and constants
sourceFile = ...
    '../data/waterdata2lambda_2Bubbles_z.mat';
data = load(sourceFile);
% trgt=data.trgt;
trgt=2;
dx=data.dx*1e3;%elements distance in [mm]
% c = mean(data.domain.ctot(data.excit_loc(1,1),data.excit_loc(1,2):end));% the images is the same for all simulations
c =1480;% % Speed of Sound [m/s]
locx = data.excit_loc(trgt,3,1)*dx;%position of the porb in mm
locy = data.excit_loc(trgt,1,1)*dx;%position of the porb in mm
locz = dx*(data.yDim-5)-data.excit_loc(trgt,2,1)*dx;%position of the porb in mm
dt =data.t(2);

p=(squeeze(data.aedata));
data = rmfield(data,'aedata'); % remove p data from structure to save memory
ss=size(p);
ychannels=ss(1);
xchannels=ss(2);
zchannels=256;
xtr=dx*(xchannels*floor(ss(2)/xchannels));%in mm
ytr=dx*(ychannels*floor(ss(1)/ychannels));%in mm

rcn1=1;%crop the time data

rcn_x2=size(p,2);%
lft_el=-xtr/2;
rght_el=xtr/2;

rcn_y2=size(p,1);
dwn_el=-ytr/2;
up_el=ytr/2;

near_coor=0;
far_coor=dx*data.yDim;

% p=p(:,:,rcn1:end);
t =data.t(rcn1+length(data.t)-ss(3):end);%time in seconds t(end)*c0*1000=distance in mm

% figure (10)
% [x,y,z] = meshgrid(1:data.zDim,1:data.xDim,1:length(t));
% v = p;
% xslice = [1,double(data.excit_loc(2,3))]; yslice =[double(data.excit_loc(2,1)),data.xDim]; zslice =[1,length(t)];
% h=slice(1e3*dx*x,1e3*dx*y,1e6*dt*z,v,round(1e3*dx*xslice),1e3*dx*yslice,1e6*dt*zslice);
% set(h,'edgecolor','none','FaceAlpha',0.75)
% axis([1e3*dx 1e3*dx*data.zDim 1e3*dx 1e3*dx*data.xDim 1e6*dt*rcn1 1e6*dt*(length(data.t))])
% view([85 15]);
% colormap jet(100)
% xlabel('Distance (mm)','FontWeight','bold')
% ylabel('Distance (mm)','FontWeight','bold')
% zlabel('Time (usec)','FontWeight','bold')
% title('Passive Emissions Detector','Color',[0 0 0],'FontSize', 22); drawnow;
%% 3D ASA
xlength=1e-3*xtr;%pixel distance this is very important
ylength=1e-3*ytr;%pixel distance this is very important
zlength=1e-3*far_coor;%image dapth
fc=0.8e6;
xx=linspace(-2*xlength/4,2*xlength/4,xchannels);
yy=linspace(-2*ylength/4,2*ylength/4,ychannels);
w=2*pi*fc;
% p=p;%p data for 128 channels
zz=linspace(0,zlength,zchannels);% this can be of any size not necessarilly equal to channels
dz=zlength/zchannels;
% find the position of fc
fk =(1/(dt))*linspace(0,1,length(t));
[gg, ff]=min(abs(fk-fc));% ff the data at fc
% extract the complex data P from each channel at fc
tic
aa=fft(p,[],3);

% convert data to k-space and shift
P=fftshift(fft2(aa(:,:,ff))); %convert P into K-space

% image reconstruction algorithm
kx = repmat(2*pi*(-floor(rcn_x2/2):ceil(rcn_x2/2)-1)/(xx(2)-xx(1))/rcn_x2,rcn_y2,1); % vector version 197 649
ky = repmat(2*pi*(-floor(rcn_y2/2):ceil(rcn_y2/2)-1)/(yy(2)-yy(1))/rcn_y2,rcn_x2,1)'; % vector version 197   649
Pv = repmat(P,1,1,zchannels);  % P matrix replicated over z: 197x649x256
t1 = sqrt(w.^2./c.^2-kx.^2-ky.^2); % this is 197   649  
asa = permute(ifft2(permute(ifftshift(Pv.*exp(1i.*reshape(t1(:)*zz,[rcn_y2 rcn_x2 zchannels])),2),[2 1 3])),[3 2 1]);
asa=2*(abs(asa)).^2;
% asa=circshift(asa, [0,0, (xchannels+1)/2]);
toc

figure (2)
subplot(3,3,[1 2 4 5])
imagesc((xx-xx(1))*1e3-locx,zz*1e3-locz,squeeze(asa(:,data.excit_loc(trgt,1,1),:)))
ylabel('Axial (mm)','FontWeight','bold')
axis image
caxis([0 0.03])

% hb = colorbar;
subplot(3,3,[3 6])
imagesc((yy-yy(1))*1e3-locy,zz*1e3-locz,squeeze(asa(:,:,data.excit_loc(trgt,3,1))))
xlabel('Elevation (mm)','FontWeight','bold')
axis image
caxis([0 0.03])

% hb = colorbar;
subplot(3,3,[7 8])
imagesc((xx-xx(1))*1e3-locx,(yy-yy(1))*1e3-locy,squeeze(asa(round((locz)/(zz(2)*1e3)),:,:)))
xlabel('Transverse (mm)','FontWeight','bold')
ylabel('Elevation (mm)','FontWeight','bold')
colormap jet;
axis image
caxis([0 0.03])


figure (4)
[x,y,z] = meshgrid(1:ychannels,1:zchannels,1:xchannels);
v = p;
xslice = (data.excit_loc(trgt,1,1)); yslice =locz; zslice =(data.excit_loc(trgt,3,1));
h=slice(dx*x,1e3*dz*y,dx*z,asa,dx*round(xslice),yslice,dx*zslice);
set(h,'edgecolor','none','FaceAlpha',0.50)
view([-80 25]);
xlabel('Distance (mm)','FontWeight','bold')
ylabel('Distance (mm)','FontWeight','bold')
zlabel('Distance (mm)','FontWeight','bold')
title('3D ASA','FontWeight','bold')