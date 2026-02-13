function Sfb = correlationMatFB(rx, L)
    % estimate the correlation matrix using forward-backward averaging 
    % similar to spatial smoothing for angle of arrival estimation
    % originally we can obtain the
    N = length(rx);
    M = N - L + 1;
    J = obtainExchangeMatrix(M);
    Sfb = zeros(M);
    for i = 1:L
        Xi = rx(i:i+M-1);
        Sfb = Sfb + (Xi'*Xi + J*transpose(Xi)*conj(Xi)*J);
    end
    Sfb = Sfb/(2*L);
end