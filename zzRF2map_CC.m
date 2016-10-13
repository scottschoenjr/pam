function [img_frm DC_bias ray] = zzRF2map_CC(rf,vox,Temperature,c)
% CC changes: nupsample=1, passed ray as output to use element positions
% elsewhere

Nupsample = 1; % normal sampling f is 4 times the central frequency (Nupsample=1)
% temp
ray = zzray('64');
ray.fs = ray.fs*Nupsample;
% ray.c=SpeedOfSound(Temperature)/1e3;
ray.c=c/1e3;
rf = zzupsample(rf,Nupsample);% upsamples the data by "Nupsample" times...makes data smother whne finite resolution is used

% variables that will be saved
img_frm = zeros(vox.Nz,vox.Nx);
DC_bias = zeros(vox.Nz,vox.Nx);

%% generate img and res
fprintf('   ');
zmax = max(vox.z);
rf_ = rf';

for Ivox = 1:vox.N
    if mod(Ivox,10)==0
        fprintf('\b\b\b%3.0f',100*Ivox/vox.N);
    end
    
    dels = zzpos2del(vox.r(Ivox,:),ray);%convert from distance of the elements (64) to image cordinates to delays
    delrf = zzdelRF(rf_,dels,Nupsample);%  rf data with different time delays
    
    %     sqdels = sqrt(dels);%gain function (because of elevational focussing, assume circular spreading)
    sqdels = ones(length(dels),1)';% no gain
    img_frm(Ivox) = sum((sqdels*delrf).^2 - (sqdels*delrf.^2));
    %     if length((sqdels*delrf).^2 - (sqdels*delrf.^2))<2.0
    %         img_frm(Ivox) = sum((sqdels*delrf).^2 - (sqdels*delrf.^2));
    %     else
    %         img_frm(Ivox) = max((sqdels*delrf).^2 - (sqdels*delrf.^2));
    %     end
    %    img_frm(Ivox) = var((sqdels*delrf).^2 - (sqdels*delrf.^2));
    %    img_frm(Ivox) = sum((dels*delrf).^2 - (dels*delrf.^2));
    %    DC_bias(Ivox) = sum((sqdels*delrf.^2));
    %    DC_bias(Ivox) = sum((ones(1,64)*delrf.^2));
    %    img_frm(Ivox) = sum((sqdels*delrf).^2);
end
fprintf('*');