function [psn, psm, nOut] = music1D(rx, SCS, fc, musicParams)
    % estimate the 1D pseudospectrum from a 2D signal in both dimensions
    % rx            :       2D signal
    % SCS           :       subcarrier spacing
    % fc            :       carrier frequency
    % musicParams   :       structure containing the configuration of the
    %                       music method

    N = musicParams.N;
    M = musicParams.M;
    Ln = musicParams.Ln;
    Lm = musicParams.Lm;
    n = musicParams.Order;
    dSearchSpace = musicParams.dSearchSpace;
    vSearchSpace = musicParams.vSearchSpace;


    rxN = transpose(rx);
    if musicParams.corrMatMethod == 0
        Sxn = zeros(N);
        Sxm = zeros(M);
        for i = 1:M
            Sxn = Sxn + rxN(i, :)'*rxN(i, :);
        end        
        for i = 1:N
            Sxm = Sxm + rx(i, :)'*rx(i, :);
        end
    elseif musicParams.corrMatMethod == 1
        Sxn = zeros(N-Ln+1);
        Sxm = zeros(M-Lm+1);
        for i = 1:M
           Sxn = Sxn + correlationMatFB(rxN(i, :), Ln); 
        end
        
        for i = 1:N
            Sxm = Sxm + correlationMatFB(rx(i, :), Lm);
        end
    else
        error('Wrong method selected')
    end


    [Un, Pn, ~] = svd(Sxn);
    nEigsN = detectEigenvalueJump(diag(Pn), 100);
    nN = min([n, nEigsN]);
    Un = Un(:, nN+1:end);
    [Um , Pm, ~] = svd(Sxm);
    nEigsM = detectEigenvalueJump(diag(Pm), 100);
    nM = min([n, nEigsM]);
    Um = Um(:, nM+1:end);

    [N, ~] = size(Sxn);
    [M, ~] = size(Sxm);

    % obtain pseudospectrum for the rows to estimate distance
    psn = zeros(1, length(dSearchSpace));
    for i = 1:length(dSearchSpace)
        baseDistance = exp(1j*2*pi*2*dSearchSpace(i)/physconst('LightSpeed')*SCS*(0:N-1)');
        psn(i) = abs(1/(baseDistance'*(Un*Un')*baseDistance));
    end

    % obtain pseudospectrum for the columns to estimate velocity
    psm = zeros(1, length(vSearchSpace));
    for i = 1:length(vSearchSpace)
        baseVelocity = exp(1j*2*pi*(-2*vSearchSpace(i)*fc/physconst('LightSpeed'))/SCS*(0:M-1)');
        psm(i) = abs(1/(baseVelocity'*(Um*Um')*baseVelocity));
    end
    nOut = max(nN, nM);
end

function ix = detectEigenvalueJump(eigs, threshold)
    % detects a jump in the values of the array eigs when the ratio between
    % two consecutive elements is higher than threshold
    % Input Parameters:
    % eigs          :           array to search for jump
    % threshold     :           threshold for the ratio between values
    ix = 1;
    while ix < length(eigs) & eigs(ix)/eigs(ix+1) < threshold
        ix = ix + 1;
    end
end