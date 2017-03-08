% plot( 1E3.*(0.1-z), fliplr(asamap( :, middleIndex )./max(abs(asamap(:, middleIndex ) ) )), 'k' );
% xlim( [20, 100] );
% xlabel( 'Axial Distance [mm]' );
% ylabel( 'Normalized Amplitude' );
% 
% figure(); 
% pcolor(zAsaPlot.*1E3, xAsaPlot.*1E3,mod(real(phi), 2.*pi ) ); 
% shading flat;
% xlabel('Axial Distance [mm]');
% ylabel('Transverse Distance [mm]');
% cbh = colorbar;
% caxis( [0, 2.*pi] );
% set( cbh, 'YTick', [0, pi./2, pi, 3.*pi./2, 2.*pi], ...
%     'YTickLabel', {'$0$';'$\pi/2$';'$\pi$';'$3\pi/2$';'$2\pi$'}, ...
%     'TickLabelInterpreter', 'LaTeX', 'FontSize', 18);
% ylim( [-40, 40] );

figure(); 
plot(z, mu_z.^(2).*c0, 'k' ); 
shading flat;
xlabel('Axial Distance [mm]');
ylabel('Sound Speed [m/s]');
ylim( [1000, 3000] );