function ray = zzray(option)
% function ray = zzray(option)
%
% return x values for Zonare array given some options
% use mm

% the wrapping 'xs_full(Ilmn([33:end 1:32]))'
% has been proven to be valid experimentally

if ~exist('option','var')
    help zzray;
    disp('options are: 64,128,FULL');
    return;
end

dx =  2*0.6412; % spacing between ACU382 elements (mm)
% positions of all the elements
% xs_full = (-dx*31.5):dx:(dx*31.5);

switch option
    case '64'
        xs_ray = [((-dx*31.5):dx:0),(0:dx:(dx*31.5))];
    case '128'
        xs_ray = (-dx*63.5):dx:(dx*63.5);
    otherwise
        xs_ray = (-dx*31.5):dx:(dx*31.5);
        
end
ray = x2ray(xs_ray);
