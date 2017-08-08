% Plot axial profile
clear all
close all
clc;

% Specify data files
dataFiles = {'layers.mat'; 'stratified.mat'; 'stratifiedSubsection.mat' };

colors = { [116, 0, 83]./255; ... % Whistle (purple)
    [249, 94, 16]./255; ... % Horizon (orange)
    [249, 94, 16]./255 ...
    };

% Plot each
figure();
hold on;
for fileCount = 1:length( dataFiles );
    load( dataFiles{ fileCount } );
    if fileCount == 3
        plot( 1E3.*z, axialProfileNorm, '--', ...
            'Color', colors{fileCount}, ...
            'LineWidth', 4 );
    else
        plot( 1E3.*z, axialProfileNorm, ...
            'Color', colors{fileCount}, ...
            'LineWidth', 3 );
    end
end
set( gca, 'XDir', 'Reverse' );
xlabel( 'Distance from Receiver [mm]', 'FontSize', 26 );
ylabel( 'Normalized Profile', 'FontSize', 26 );
box on;


