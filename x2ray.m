%**************************************************************************
%
% Convert x to Ray
%
%   Function converts xs vector into a ray struct.
%
% Inputs
%   xs - Array of x-positions
%
% Outputs
%   ray. - Struct with fields
%     fs - Sampling frequency [Hz]
%     N  - Number of beams
%     z  - z-coordinate of beam [mm]
%     y  - y-coordinate of beam [mm]
%     r  - Ray vector (3-by-1) [mm]
%     A  - Magnitude of vector [m]
%
%**************************************************************************

function ray = x2ray(xs)

% function vox = x2ray(xs)
% converts vector of 'xs' in mm into 'ray' structure
% NEW STANDARD UNITS

% convert into m
% xs = xs/1000;

% Make sure z runs along col
if size(xs,1) == 1
    xs = xs';
end

% Set sampling frequency
ray.fs = 3.21*4; % [MHz] Verasonics for RF (MHz) make it general
ray.x = xs;
ray.N = length(ray.x);
ray.z = zeros(ray.N,1);
ray.y = zeros(ray.N,1);

ray.r = [ray.x ray.y ray.z];
ray.A = norm(xs(end)-xs(1));

end