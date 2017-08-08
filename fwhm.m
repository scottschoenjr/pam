% *************************************************************************
% Find Maximum and FWHM of 1D data trace
%
%   fwhm calculates the Full Width at Half Maximum (FWHM) of a
%   positive 1D input function yData(xData) with spacing given by xData.
%
% Inputs
%   yData    - 1-D data data
%   yData    - xData data (or dx)
%   plotFwhm - [optional 1 to plot result or 0 not to
%
% Outputs
%   fwhm      - FWHM of data
%   xMaxValue - Coordinate of maximum value
%
% Based on fwhm.m from k-Wave Toolbox (http://www.k-wave.org),
% © 2009-2014 Bradley Treeby and Ben Cox, licensed under GNU/GPL v3.
% 
%**************************************************************************

function fwhm = fwhm(yData, xData, plotFwhm)

% check for optional inputs
if nargin == 2
    plotFwhm = 1;
end

% check if dx is given in place of an xData array
if length(xData) == 1
    xData = 0:xData:(length(yData)-1)*xData;
end

% check the input is 1D
if numDim(yData) ~= 1
    error('input function must be 1 dimensional');
end

% find the maximum value of yData(xData)
[f_max, i_max] = max(yData);

% setup the indexing variables for finding the leading edge
index = i_max;

% loop until the index at half maximum is found
while yData(index) > 0.5*f_max 
    index = index - 1;
    if index < 1
        error('left half maximum not found');
    end
end

% linearly interpolate between the previous values to find the value of xData
% at the leading edge at half maximum
m = (yData(index+1) - yData(index))/(xData(index+1) - xData(index));
c = yData(index) - xData(index)*m;
i_leading = (0.5*f_max  - c)/m;

% setup the indexing variables for finding the trailing edge
index = i_max;

% loop until the index at half maximum is found
while yData(index) > 0.5*f_max 
    index = index + 1;
    if index > length(yData)
        error('right half maximum not found');
    end
end

% Linearly interpolate between the previous values to find the value of 
% xData at the trailing edge at half maximum
m = (yData(index-1) - yData(index))/(xData(index-1) - xData(index));
c = yData(index) - xData(index)*m;
i_trailing = (0.5*f_max  - c)/m;
    
% compute the FWHM
fwhm = abs(i_trailing - i_leading);

% plot the function and the FWHM if required
if plotFwhm   
    figure;
    plot(xData, yData, 'b-');
    xlabel('x');
    ylabel('f(x)');    
    hold on;
    plot([i_leading i_trailing], [0.5*f_max 0.5*f_max], 'r*-');
    title(['FWHM: ' num2str(fwhm), ' at x = ', num2str(xData(i_max)) ]);
end

end