%**************************************************************************
%
% Get Delays
%
%   Function returns the delays between for the input position and ray
%   object.
%
% Inputs
%   pos
%   ray
%
% Output
%   dels
%
%
%**************************************************************************


function dels = zzpos2del(pos,ray)

dels = round( ...
    sqrt(sum((ones(64,1) * pos -ray.r).^2,2)) * ray.fs/ray.c)';
end

