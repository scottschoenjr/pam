function vox = xz2vox(xs,zs)

% function vox = xz2vox(xs,zs)
% converts vector of 'x', 'z' in mm into 'vox' structure

[X Z] = meshgrid(xs,zs);

vox.x = xs;
vox.z = zs;

vox.Nx = length(vox.x);
vox.Nz = length(vox.z);

vox.N = vox.Nz*vox.Nx;
vox.r = zeros(vox.N,3);

for Ivox = 1:vox.N
Ix = ceil(Ivox/vox.Nz);
Iz = mod(Ivox-1,vox.Nz)+1;
    vox.r(Ivox,:) = [vox.x(Ix) 0 vox.z(Iz)]; % voxels go down columns
end