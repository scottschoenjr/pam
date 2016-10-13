%**************************************************************************
%
% Get Distance to Voxel
%
%   Function returns the distance between the input position and the 
%
%measures the distances of the elements from the image coordinates and
% converts them to wavelength in order to determine the number of delays 


function dels = zzpos2del(pos,ray)

dels = round( ...
    sqrt(sum((ones(64,1) * pos -ray.r).^2,2)) * ray.fs/ray.c)';

