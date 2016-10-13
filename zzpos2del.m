function dels = zzpos2del(pos,ray)
% measures the distances of the elements from the image coordinates and
% converts them to wavelength in order to determine the number of delays 
dels = round( ...
    sqrt(sum((ones(64,1) * pos -ray.r).^2,2)) * ray.fs/ray.c)';

