function res = zzupsample(in,factor,Bcol)

% function res = zzupsample(in,factor,Bcol)
%
% use interp to upsample each row/column (column iff Bcol)
%
% Bcol=0 if not specified

if ~exist('Bcol','var')
    Bcol = 0;
end


res = zeros(size(in).*[factor 1]);

if Bcol
    res = zeros(size(in).*[factor 1]);
    for Icol=1:size(in,2)
        res(:,Icol) = interp(in(:,Icol),factor);
    end
else
    res = zeros(size(in).*[1 factor]);
    for Irow=1:size(in,1)
        res(Irow,:) = interp(in(Irow,:),factor);    
    end
end