%% from Costas Arvanitis (C) August 2010 BWH
%passive ultra-sonography algorithm based on Back-propagation
%it reconstructs cavitation selective images and superimposes them in order to constract maps of cavitation activity .
%The spectral response of the cavitation events (temporal inforamtion) and of the cavitation map (treatment infromation)are also extracted.
close all;
clear all;
%% parameters and constants
skul_cor='y'
tic
Sdir='C:\Users\carvanitis7\OneDrive - Georgia Institute of Technology\FDTD-ASA-Scott\data\simulated\Human';
Sdir1='880kHz_aeattscabs_2Dhead_crop_res100um_1gps';
data=load([Sdir filesep ([Sdir1 '.mat'])]);
% trgt=4;
dx=data.dx*1e3;%elements distance in [mm]
locz = dx*(data.xDim-5)-data.excit_loc(1,2)*dx;%position of the porb in mm
% locx=(data.excit_loc(trgt,1)+2*round(data.ratio*data.num))*dx;
locx=(data.excit_loc(1,1))*dx;
distance=round(data.ratio*data.num)*dx;
dt =data.t(2);
t = data.t;%time in seconds t(end)*c0*1000=distance in mm

rfdata=data.aedata;
ss=size(rfdata);
channels=64;
xtr=dx*(channels*floor(ss(1)/channels));%in mm
rcn_x1=3500;
rcn_x2=size(rfdata,2);
lft_el=-xtr/2;
rght_el=xtr/2;
step=0.75;%in images is 0.5;
near_coor=60;
% far_coor=dx*data.xDim;
far_coor=169;
parameters=[channels;xtr;rcn_x1;rcn_x2;lft_el;rght_el;step;near_coor;far_coor];
rf128=rfdata(:,rcn_x1:rcn_x2);
vox = xz2vox(lft_el:step:rght_el,near_coor:step:far_coor);% image coordinates of the reconstracted image
shift=dx*(size(rf128,1)-floor(size(rf128,1)/channels)*channels-1);%this is to correct for the location of the point source after binning
rf=zeros(channels,rcn_x2-rcn_x1+1);
for jj=1:channels
    rf(jj,:)=sum(rf128((floor(ss(1)/channels))*(jj-1)+1:(floor(ss(1)/channels))*jj,:));
end
rf=rf./sqrt(floor(ss(1)/channels));%correct for different binning
%display binned rf data
figure (2)
imagesc(1e6*t,0:channels-1,rf)
ylabel('Channel','FontWeight','bold')
xlabel('Time (usec)','FontWeight','bold')
hb = colorbar;
ylabel(hb,'Pressure (Pa)','FontSize',12,'FontWeight','bold')

%% passive ultrasonography algorithm
%         ttToptions%for quantitative imaging
%         im= ttRF2map(rf(:,rcn_x1:rcn_x2),options);%for quantitative imaging
vox = xz2vox((lft_el:step:rght_el),near_coor:step:far_coor);% for zooming at specific point in the reconstracted image

%% Speed of sound
if skul_cor=='y'
    ctsos=flipdim(data.domain.ctot(5:end,:),2);
    mapsos=zeros(channels,floor((data.domain.dx/step)*(data.xDim-5)));
    [ch01 dpth]=size(ctsos(1:(floor(data.yDim/channels)),1:floor(step/data.domain.dx)));
    for jj=1:channels
        for ii=1:floor((data.domain.dx/step)*(data.xDim-5))
            mapsos(jj,ii)=sum(sum(ctsos((floor(data.yDim/channels))*(jj-1)+1:(floor(data.yDim/channels))*jj,...
                (ii-1)*dpth+1:ii*dpth)))/(ch01*dpth);
        end
    end
    sos=mapsos;
    %             mapsos(:,end-12/step:end)=1479;
    %             sos=[mapsos(:,end-12/step:end),mapsos(:,1:end-12/step+1)];
else
    sos = mean(data.domain.ctot(data.excit_loc(1,1),data.excit_loc(1,2):end-5));% the images is the same for all simulations
%                 sos =1479.0;% % Speed of Sound [m/s]
end
%% passive ultrasonography algorithm
del = xtr/channels; % for xx channles
xs = [((-del*channels/2):del:0),(del:del:(del*(channels/2-1)))];
ray.fs = 1e-6/dt; % sampling rate (MHz)
ray.x = xs';
ray.N = length(ray.x);
ray.z = zeros(ray.N,1);
ray.y = zeros(ray.N,1);
ray.r = [ray.x ray.y ray.z];
% variables that will be saved
img_frm = zeros(vox.Nz,vox.Nx);
DC_bias = zeros(vox.Nz,vox.Nx);
%% generate img and res
fprintf('   ');
zmax = max(vox.z);
rf_ = rf';
avsos=zeros(1,channels);
for Ivox = 1:vox.N
    if mod(Ivox,10)==0
        fprintf('\b\b\b%3.0f',100*Ivox/vox.N);
    end
    %% sos
    if skul_cor=='y'
        [val xxm]=min(abs(ones(channels,1) * vox.r(Ivox,1)-ray.r(:,1)));
        yym=vox.r(Ivox,3)/step+1;
        xi=[ones(channels,1) ones(channels,1)*yym];
        yi=[(1:channels)' ones(channels,1)*xxm];
        for kk=1:channels
            avsos(1,kk) = mean(improfile(sos,xi(kk,:),yi(kk,:)));
        end
        dels = round(sqrt(sum((ones(channels,1) * vox.r(Ivox,:) -ray.r).^2,2)) *1000*ray.fs./avsos')';
    else
        dels = round(sqrt(sum((ones(channels,1) * vox.r(Ivox,:) -ray.r).^2,2)) *1000*ray.fs/sos)';
    end
    [NRF Nchan] = size(rf'); % in fact this is rf'
    
    Ddel = max(dels)-min(dels); % maximum difference in delays
    dels = dels-min(dels) + (0:(Nchan-1))*NRF;
    Nsample=1;
    Ip = 1:Nsample:(NRF-Ddel);
    NIp = length(Ip);
    
    delrf = rf_(ones(Nchan,1)*Ip + dels' * ones(1,NIp) );
    sqdels = sqrt(dels);% gain function (because of elevational focussing, assume circular spreading)
    %             sqdels = ones(length(dels),1)';% no gain
    if length((sqdels*delrf).^2 - (sqdels*delrf.^2))<2.0
        img_frm(Ivox) = sum((sqdels*delrf).^2 - (sqdels*delrf.^2));
    else
        img_frm(Ivox) = max((sqdels*delrf).^2 - (sqdels*delrf.^2));
    end
end
fprintf('*');
toc
figure
subplot(1,2,1)
imagesc(vox.x,vox.z,img_frm)
axis('image');
colormap jet;
hb = colorbar;
ylabel(hb,'a.u','FontSize',16,'FontWeight','bold')
ylabel('Transverse (mm)','FontWeight','bold')
xlabel('Axial (mm)','FontWeight','bold')

%align image to simulation parameters

pus=imrotate(img_frm,-90);
% fpus2=flipdim(pus,1);
uszz1=vox.z;
% usxx=vox.x-min(vox.x);
% uszz=vox.z-shift;
usxx1=vox.x-min(vox.x);
Nup=4;
usxx = interp(usxx1,Nup);
uszz = interp(uszz1,Nup);
fpus = interp2(uszz1,usxx1',pus,uszz,usxx');
subplot(1,2,2)
imagesc((locz-uszz),locx-usxx,flipdim(fpus,2))
set(gca,'YDir','Normal')
set(gca,'XDir','Normal')

axis('image');
colormap jet;
hb = colorbar;
ylabel(hb,'a.u','FontSize',16,'FontWeight','bold')
ylabel('Axial (mm)','FontWeight','bold')
xlabel('Transeverse (mm)','FontWeight','bold')

figure (102)
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
fig3 = figure(102);
position = get(fig3,'Position');
outerpos = get(fig3,'OuterPosition');
borders = outerpos - position;
edge = borders(1);
pos1 = [edge,scnsize(4) * (1/10),scnsize(3)/0.9 + edge,scnsize(4)/2];
set(fig3,'OuterPosition',pos1);
subplot(1,3,1)
imagesc((locz-uszz),locx-usxx,flipdim(fpus,2));
set(gca,'YDir','Normal')
set(gca,'XDir','Normal')
if skul_cor=='y'
    title('Passive Mapping (Variable SOS)','FontSize',12,'FontWeight','bold')
else
    title('Passive Mapping (Average SOS)','FontSize',12,'FontWeight','bold')
end
axis('image');
colormap jet;
hb = colorbar;
ylabel('Transverse (mm)','FontSize',12,'FontWeight','bold')
xlabel('Axial (mm)','FontSize',12,'FontWeight','bold')

subplot(1,3,2)
[xx1 yy1]=max(max(fpus'));
plot (flipdim((locz-uszz),2),fpus(yy1,:)/xx1,'k','LineWidth',3);title('Axial Line Profile','FontSize',12,'FontWeight','bold')
% xlim([-40,40])
ylim([0,1.0])
xlabel('Distance From Source (mm)','FontSize',12,'FontWeight','bold')
ylabel('Normalized to Max','FontSize',12,'FontWeight','bold')
% line([0,0],[0, 1.2],'color',[0 0 0],'LineWidth',2);
% find the fwhm
aa=sort(abs(fpus(yy1,:)/xx1-0.5));
zz=flipdim((locz-uszz),2);
zznr=zz(abs(fpus(yy1,:)/xx1-0.5)==aa(1));
zzfr=zz(abs(fpus(yy1,:)/xx1-0.5)==aa(2));
line([zznr,zznr],[0, 1.2],'color',[1 0 0],'LineWidth',2);
line([zzfr,zzfr],[0, 1.2],'color',[1 0 0],'LineWidth',2);
text(ceil(abs(zznr))+1, 0.95,['FWHM ' num2str(abs(zznr-zzfr)) ' mm'],'Rotation',0,'FontSize',12,'FontWeight','normal','Color','k')

subplot(1,3,3)
[xx2 yy2]=max(max(fpus));
plot (locx-usxx,fpus(:,yy2)/xx2,'k','LineWidth',3);title('Transverse Line Profile','FontSize',12,'FontWeight','bold')
xlim([-8,8])
ylim([0,1.0])
xlabel('Distance From Source (mm)','FontSize',12,'FontWeight','bold')
ylabel('Normalised to Max','FontSize',12,'FontWeight','bold')
% line([0,0],[0, 1.2],'color',[0 0 0],'LineWidth',2);
aa=sort(abs(fpus(:,yy2)/xx2-0.5));
xx=locx-usxx;
xxnr=xx(abs(fpus(:,yy2)/xx2-0.5)==aa(1));
xxfr=xx(abs(fpus(:,yy2)/xx2-0.5)==aa(2));
line([xxnr,xxnr],[0, 1.2],'color',[1 0 0],'LineWidth',2);
line([xxfr,xxfr],[0, 1.2],'color',[1 0 0],'LineWidth',2);
text((abs(xxfr)), 0.95,['FWHM ' num2str(abs(xxfr-xxnr)) ' mm'],'Rotation',0,'FontSize',12,'FontWeight','normal','Color','k')

fwhm=[abs(zznr-zzfr); abs(xxfr-xxnr)];

% %% Save data
% aa=regexp(Sdir,'\');
% folder1=[Sdir,'\analysis'] ;
% aa1=regexp(Sdir1,'_');
% folder2=['PAM_ncrop' Sdir1(1:aa1(1)-1) Sdir1(aa1(2):aa1(3)-1) '_abs_res' num2str(1000*step) '_skul_cor_' (skul_cor) Sdir1(aa1(end):end) '_02'];
% savefile =[folder1 '\' folder2 '.mat'];
% save(savefile,'fpus','uszz','usxx','locz','locx','shift', 'fwhm','parameters');
% saveas(102,[folder1 '\' folder2], 'fig');