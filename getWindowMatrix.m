function W = getWindowMatrix(N, M, window)
    if strcmp(window, 'hamming')
        wN = hamming(N);
        wM = hamming(M)';
    elseif strcmp(window, 'blackman')
        wN = blackman(N);
        wM = blackman(M)';
    elseif strcmp(window, 'chebyshev')
        wN = chebwin(N);
        wM = chebwin(M)';
    elseif strcmp(window, 'rect')
        wN = ones(N, 1);
        wM = ones(1, M);
    else
        error('Wrong window type');
    end
    % W = 1/(norm(wN)*norm(wM))*wN*wM;
    W = wN*wM;
end