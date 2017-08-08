% Script to plot phase corrections for stratified ASA 

figure();
plot(1E3.*z, c_z./c, 'k');
set( gca, 'XDir', 'reverse' );
xlabel( 'Distance from Receiver [mm]', 'FontSize', 26 );
ylabel( '$$c/c_{\rm avg}$$', 'FontSize', 28 );
ylim( [0.5, 1.5] );
xlim( [0, 80] );

figure();
subplot( 4, 1, 4)
plot(1E3.*z, c_z./c, 'k');
set( gca, 'XDir', 'reverse' );
xlabel( 'Distance from Receiver [mm]', 'FontSize', 26 );
ylabel( '$$c/c_{\rm avg}$$', 'FontSize', 24 );
ylim( [0.5, 1.5] );
xlim( [0, 80] );

subplot( 4, 1, 1:3)
pcolor( zAsaPlot.*1E3, xAsaPlot.*1E3, mod( phi, 2.*pi ) );
set( gca, 'XDir', 'reverse', 'XTIckLabel', '' );
ylabel( 'Transverse Distance [mm]', 'FontSize', 22 );
ylim( [-40, 40] );
shading interp;