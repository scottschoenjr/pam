function [rim rvox] = zzUP(im,Nup,vox)


if ~exist('vox','var')
    vox = xz2vox(1:size(im,2),1:size(im,1));
end

if isscalar(Nup)
    rx = interp(vox.x,Nup);
    rz = interp(vox.z,Nup);
    rvox = xz2vox(rx,rz);

    rim = interp2(vox.x,vox.z',im,rx,rz');
end

if nargout==0
%     imagesc(rx,rz,nm(rim));
    imagesc(rx,rz,(rim));%not normalised
%     colorbar
end
