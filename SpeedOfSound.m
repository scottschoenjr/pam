%% speed of sound calculation

function sos=SpeedOfSound(Temperature)

sos=1402.385 + 5.038813*Temperature - 0.05799136 *(Temperature)^2 ...
+ 0.0003287156*(Temperature)^3 - 0.000001398845*(Temperature)^4 ...
+ 0.00000000278786*(Temperature)^5;