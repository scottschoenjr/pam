%**************************************************************************
%
% Add White Gaussian Noise
%
%   Function adds white Gaussian noise to the input signal with the input
%   SNR.
%
% Inputs
%   signal - Vector of values to which to add noise. If signal is a matrix,
%            noise is added along the rows. If the values are complex, then
%            noise is added separately to the real and imaginary parts.
%   SNR_dB - Signal to noise ratio desired. If a vector, signal must be a
%            matrix with the same number of rows.
%
% Outputs
%   noisySignal - Noisy signal.
%   result      - 0 if successful and an error string otherwise
%                 (noisySignal will be NaN)
%
% References
%    Consulted Mathuranathan Viswanathan's "add_awgn_noise()".
%
%**************************************************************************

function [ noisySignal, result ] = addGaussianWhiteNoise( signal, SNR_dB )

% Initialize
noisySignal = NaN;
result = '';

% Check inputs

% Make sure of maximum dimension 2
signalDimensions = size(signal);
if length( signalDimensions ) >= 3
    result = 'Signal must be 1- or 2-D';
    return;
end
 
% Determine if 1-D and make sure it has the correct dimension (and flag so
% we can switch it back)
signalIs1D = any( signalDimensions == 1);
inputAsColumnVector = 0;
if signalIs1D
    if signalDimensions(1) ~= 1
        % Make row vector
        inputAsColumnVector = 1;
        signal = signal'; 
    end
end
numRows = size( signal, 1 );
numCols = size( signal, 2 );

% Make sure the correct number of SNRs is specified
numSNRs = length(SNR_dB);
validNumSNRs = ...
       ( numSNRs == 1 ) ...
    || ( signalDimensions(1) == numSNRs );
if ~validNumSNRs
    result = 'Must specify one SNR or one SNR for each row.';
    return;
end

% Add in noise to each row
noisySignal = 0.*signal;
for rowCount = 1:numRows
    
    % Get current signal and SNR
    currentSignal = signal( rowCount, : );
    if numSNRs > 1
        SNR = 10.^(SNR_dB(rowCount)./10); % Linear SNR
    else
        SNR = 10.^(SNR_dB./10); % Linear SNR
    end
    totalEnergy = sum(abs(currentSignal).^(2))./length(currentSignal);
    spectralDensity = totalEnergy./SNR;
    
    % Determine if signal is complex, and separate into real and imaginary
    % if it is parts
    signalComplex = any( ~isreal(currentSignal) );
    if signalComplex
        noiseSigma = sqrt( spectralDensity./2 );
        noise = noiseSigma.*( randn(1, numCols) + 1j.*randn(1, numCols) );
    else
        % Otherwise just compute the noise
        noiseSigma = sqrt( spectralDensity );
        noise = noiseSigma.*randn(1, numCols);
    end
    
    % Add the noise into the signal
    noisySignal( rowCount, : ) = currentSignal + noise; 
    
end

% Flip back 1-D vector if requried
if inputAsColumnVector
    noisySignal = noisySignal';
end

% Return success
result = 0;

end