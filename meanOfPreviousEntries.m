%**************************************************************************
%
% Get Mean of Previous Elements
%
%   Function returns a vector whose entries are the mean of all previous
%   entries of the input vector. So, e.g., for 
%        x = [ 1, 3, 5, 3]
%   the returned vector is
%        y = [ 1, 2, 3, 3]
%   If x is a 2D array, the averaging is performed along the input 
%   dimension: 1 to average down the columns, 2 to averge along the rows.
%
% Inputs
%   x   - Array to be averaged
%   DIM - Dimension to average along [opt., default = 1]
%
% Outputs
%   xAvg - Averaged array
%
%           Scott Schoen Jr | Georgia Tech | 20170108
%
%**************************************************************************

function [ xAvg ] = meanOfPreviousEntries( x, DIM )

% Get the dimension to average along if not specified
if nargin < 2 || ~isnumeric( DIM )
    DIM = 1;
else
    DIM = round( DIM );
end

% Make sure valid dimension
numDims = length( find( size(x) > 1 ) );
if DIM > 2 % to be numDims eventually
    errorString = sprintf( ...
        'Can''t average along dimension %2d, only %2d dimensions', ...
        DIM, numDims );
    error( errorString );
    xAvg = NaN;
    return;
end
DIM = max( 1, round( DIM ) ); % In case of non-integer inputs

% We'll just average along columns, so if rows are desired, flip the aray
% here and we'll flip it back at the end.
% Flip array if we need to
flipFlag = 0;
if DIM == 2
    x = x';
    flipFlag = 1;
end

% Initialize
xAvg = 0.*x;

% Perform averaging along specified dimension

% For array
if numDims == 2
    numEntriesInDimension = length( x(:, 1) );
    xAvg( 1, : ) = x( 1, : );
    for entryCount = 2 : numEntriesInDimension
        xAvg( entryCount, : ) = mean( x( 1:entryCount, : ) );
    end
% For vectors
else
    numEntriesInDimension = length( x );
    xAvg(1) = x(1);
    for entryCount = 2 : numEntriesInDimension
        xAvg( entryCount ) = mean( x( 1:entryCount ) );
    end
end

% Flip back if necessary
if flipFlag
    xAvg = xAvg';
end

end