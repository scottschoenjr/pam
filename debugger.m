% Loop debugger
clear all
close all
clc

Nx = 100;
Ny = 100;
xx = linspace( -1, 1, Nx );
yy = linspace( -1, 1, Ny );
[x, y] = meshgrid( xx, yy );

z = x.^(2).*sin(8.*x) + cos( 7.*pi.*y ) + sin( 13.*pi.*y );

dx = xx(2) - xx(1);
dy = yy(2) - yy(1);
kkx = 2.*pi.*( -Nx/2 : Nx/2 - 1 )./(dx);
kky = 2.*pi.*( -Ny/2 : Ny/2 - 1 )./(dy);
[kx, ky] = meshgrid( kkx, kky );
Az = fftshift( fft2( z ) );

figure()
pcolor( x, y, z );
shading flat;

figure()
pcolor( kx, ky, abs(Az) );
shading flat;

figure()
pcolor( kx, ky, angle(Az) );
shading flat;


