%**************************************************************************
%
% Convert Position to Voxel
%
%   Function converts position into a voxel structure.
%
% Inputs
%   xs - x position [mm]
%   y  - y position [mm]
%
% Outputs
%   vox.  - Structure with fields
%     Nx - Number of x units
%     Ny - Number of y units
%     N  - Total number of units
%     
%
%function vox = xz2vox(xs,zs)
% converts vector of 'x', 'z' in mm into 'vox' structure

function vox = xz2vox(xs,zs)



[X Z] = meshgrid(xs,zs);


% Store input
vox.x = xs;
vox.z = zs;

% Get the number of units
vox.Nx = length(vox.x);
vox.Nz = length(vox.z);

vox.N = vox.Nz*vox.Nx;
vox.r = zeros(vox.N,3);

for Ivox = 1:vox.N
    Ix = ceil(Ivox/vox.Nz);
    Iz = mod(Ivox-1,vox.Nz)+1;
    vox.r(Ivox,:) = [vox.x(Ix) 0 vox.z(Iz)]; % voxels go down columns
end

end