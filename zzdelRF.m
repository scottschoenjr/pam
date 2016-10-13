function res = zzdelRF(rf,dels,Nsample)

% function res = zzdelRF(rf,dels,Nsample)
% delays RF lines with delays 'dels'
% Nsample is downsampling or decimation ratio (default 1)

if ~exist('Nsample','var')
    Nsample = 1;
end

[NRF Nchan] = size(rf); % in fact this is rf'

Ddel = max(dels)-min(dels); % maximum difference in delays
dels = dels-min(dels) + (0:(Nchan-1))*NRF;

Ip = 1:Nsample:(NRF-Ddel);
NIp = length(Ip);

res = rf( ones(Nchan,1)*Ip + dels' * ones(1,NIp) );
% figure (111)
% imagesc(res)