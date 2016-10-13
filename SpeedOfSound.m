%**************************************************************************
%
% Speed of Sound in Water
%
%   Function returns the speed of sound in water for the given input
%   temperature, calculated from Eq. (1) of "Speed of Sound in Pure Water",
%   De Grosso and Mader JASA (52)1442, 1972.
%
% Input
%   temp_degC - Array of water temperatures [deg C]
%
% Output
%   soundSpeed - Sound speed [m/s]
%
%**************************************************************************

function [ soundSpeed ] = SpeedOfSound( temp_degC )

% Coefficients
k0 = 1402.385;
k1 = 5.038813;
k2 = -0.05799136;
k3 = 0.0003287156;
k4 = -0.000001398845;
k5 = 0.00000000278786;

% Compute sound speed
T = temp_degC;
soundSpeed = ... % Eq.(1), values from Table III (p. 1443)
    k0 + k1.*T + k2.*T.^(2) + k3.*T.^(3) + k4.*T.^(4) + k5.*T^(5);

end