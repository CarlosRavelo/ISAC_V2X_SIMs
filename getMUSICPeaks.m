function musicPeaks = getMUSICPeaks(rx, Psn, Psm, SCS, fc, peakSearchParams)
    c = physconst('LightSpeed');

    % extract config params
    peakPromin = peakSearchParams.peakProminence;
    Ln = peakSearchParams.FBLn;
    Lm = peakSearchParams.FBLm;
    n = peakSearchParams.nOrder;
    dSpace = peakSearchParams.dSpace;
    vSpace = peakSearchParams.vSpace;

    [N, M] = size(rx);

    % search for peaks
    [peaksN, locN] = findpeaks(Psn, 'MinPeakProminence', peakPromin);
    [peaksM, locM] = findpeaks(Psm, 'MinPeakProminence', peakPromin);

    if length(locN) > 1 || length(locM) > 1
        S2dfb = correlationMat2DFB(rx, Ln, Lm);
        N2 = N - Ln + 1;
        M2 = M - Lm + 1;
        [U2, ~, ~] = svd(S2dfb);
        U2n = U2(:, n+1:end);
    
        % compare peaks on each dimmension and initialize the results
        nPeaksN = length(peaksN);
        nPeaksM = length(peaksM);
        musicPeaks = zeros(max([nPeaksN, nPeaksM]), 2);
    
        if nPeaksN > nPeaksM
            for i = 1:nPeaksN
                maxPeak = -inf;
                baseDistance = exp(1j*2*pi*2*dSpace(locN(i))/c*SCS*(0:N2-1)');
                for k = 1:nPeaksM
                    baseVelocity = exp(-1j*2*pi*(2*vSpace(locM(k))*fc/c)/SCS*(0:M2-1));
                    v2D = reshape(baseDistance*baseVelocity, N2*M2, 1);
                    Psik = abs(1/(v2D'*(U2n*U2n')*v2D));
                    if Psik > maxPeak
                        maxPeak = Psik;
                        musicPeaks(i, :) = [vSpace(locM(k)), dSpace(locN(i))];
                    end
                end
            end
        else
            for i = 1:nPeaksM
                maxPeak = -inf;
                baseVelocity = exp(-1j*2*pi*(2*vSpace(locM(i))*fc/c)/SCS*(0:M2-1));
                for k = 1:nPeaksN
                    baseDistance = exp(1j*2*pi*2*dSpace(locN(k))/c*SCS*(0:N2-1)');
                    v2D = reshape(baseDistance*baseVelocity, N2*M2, 1);
                    Psik = abs(1/(v2D'*(U2n*U2n')*v2D));
                    if Psik > maxPeak
                        maxPeak = Psik;
                        musicPeaks(i, :) = [vSpace(locM(i)), dSpace(locN(k))];
                    end
                end
            end
    
        end
    else
        musicPeaks = [vSpace(locM), dSpace(locN)];
    end
    
    
end

function S2dfb = correlationMat2DFB(rx, Ln, Lm)
    % estimate the correlation matrix of the 2D signal rx using Ln and 
    % Lm subsamples in the row and column directions respectively
    [N, M] = size(rx);
    Nl = N - Ln + 1;
    Ml = M - Lm + 1;
    J = obtainExchangeMatrix(Nl*Ml);
    S2dfb = zeros(Nl*Ml);
    for i = 1:Ln
        for k = 1:Lm
            rxI = reshape(rx(i:i+Nl-1, k:k+Ml-1), 1, Nl*Ml);
            S2dfb = S2dfb + (rxI'*rxI + J*transpose(rxI)*conj(rxI)*J);
        end
    end
    S2dfb = S2dfb/(2*i*k);
end