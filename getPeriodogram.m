function [periodogram, dAxis, vAxis] = getPeriodogram(rx, fc, SCS, varargin)
    C = physconst('LightSpeed');
    % check if the FFT and iFFT length where specified
    if nargin > 3
        periodogramParams = varargin{1};
    elseif nargin == 3
        periodogramParams = struct;
        [periodogramParams.N, periodogramParams.M] = size(rx);
        periodogramParams.EffBW = SCS*periodogramParams.N;
        periodogramParams.Window='hamming';
        periodogramParams.dMax = 70;
        periodogramParams.vMax = 30;
    else
        error('Wrong number of input arguments');
    end
    % obtain periodogram
    W = getWindowMatrix(size(rx,1), size(rx,2), ...
        periodogramParams.Window);
    P = ifft(rx.*W, 2^nextpow2(size(rx,1)), 1);
    P = fft(P, 2^nextpow2(size(rx,2)), 2);
    periodogram = abs(fftshift(P, 2).^2);
    
    distanceAxis = (0:size(P,1)-1)*C/(2*SCS*size(P,1))/periodogramParams.deltaN;
    velocityAxis = (-size(P,2)/2:size(P,2)/2 - 1)*C*SCS/(2*fc*size(P,2))/periodogramParams.deltaM;
    dIx = find(distanceAxis <= 1.5*periodogramParams.dMax);
    vIx = find(velocityAxis >= -1.5*periodogramParams.vMax & velocityAxis <= 1.5*periodogramParams.vMax);
    dAxis = distanceAxis(dIx);
    vAxis = velocityAxis(vIx);
    periodogram = periodogram(dIx, vIx);
end