function ray = x2ray(xs)

% function vox = x2ray(xs)
% converts vector of 'xs' in mm into 'ray' structure
% NEW STANDARD UNITS

% convert into m
% xs = xs/1000;

% make sure z runs along col
if size(xs,1)==1
    xs = xs';
end
ray.fs = 3.21*4; % sampling rate of Verasonics for RF (MHz) make it general
ray.x = xs;
ray.N = length(ray.x);
ray.z = zeros(ray.N,1);
ray.y = zeros(ray.N,1);

ray.r = [ray.x ray.y ray.z];
ray.A = norm(xs(end)-xs(1));