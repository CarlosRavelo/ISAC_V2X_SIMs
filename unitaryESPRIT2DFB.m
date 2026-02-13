function results = unitaryESPRIT2DFB(rx, SCS, fc, espritParams)
    c = physconst('LightSpeed');
    % extract parameters for processing
    kN = espritParams.Nlag;
    kM = espritParams.Mlag;
    nIn = espritParams.Nin;
    mIn = espritParams.Min;
    n = espritParams.Order;

    [N, M] = size(rx);
    kRows = N - kN + 1;         % different realizations on the row direction
    kColumns = M - kM + 1;      % same in column directions

    % obtain auxiliar matrices
    Qn = getQTransformationMatrix(nIn);
    Qn1 = getQTransformationMatrix(nIn - 1);

    Js1N = [eye(nIn-1), zeros(nIn-1, 1)];
    Js2N = [zeros(nIn-1, 1), eye(nIn-1)];

    K1 = real(Qn1'*Js1N*Qn);
    K2 = imag(Qn1'*Js2N*Qn);

    Qm = getQTransformationMatrix(mIn);
    Qm1 = getQTransformationMatrix(mIn-1);
    Kx1 = kron(eye(mIn), K1);
    Kx2 = kron(eye(mIn), K2);

    Js2M = [zeros(mIn-1, 1), eye(mIn-1)];

    K3 = real(Qm1'*Js2M*Qm);
    K4 = imag(Qm1'*Js2M*Qm);
    Ky1 = kron(K3, eye(nIn));
    Ky2 = kron(K4, eye(nIn));

    Y = zeros(nIn*mIn, kRows*kColumns);
    index = 1;

    nSignals = kN - nIn + 1;        % number of subsignals per realization in row
    mSignals = kM - mIn + 1;        % and column directions

    for i = 1:kRows
        for k = 1:kColumns
            rxI = rx(i:i+kN-1, k:k+kM-1);
                for ii = 1:nSignals
                    for kk = 1:mSignals
                        rxIn = rxI(ii:ii+nIn-1, kk:kk+mIn-1);
                        Yii = reshape(Qn'*rxIn*conj(Qm), nIn*mIn, 1);
                        Y(:, index) = Y(:, index) + Yii;
                    end
                end
            index = index + 1;
        end
    end

    Ymat = [real(Y), imag(Y)];

    % obtain the SVD
    [U, P, ~] = svd(Ymat);
    nEig = detectEigJumpBackwards(diag(P), 5, n);
    % if ~nEig
    %     printf('Here')
    % end
    nReal = min([n, nEig]);
    Us = U(:, 1:nReal);

    Psix = (Kx1*Us)\(Kx2*Us);
    Psiy = (Ky1*Us)\(Ky2*Us);
    PsiMatrix = Psix + 1j*Psiy;
    lambdas = eig(PsiMatrix);

    % estimate the target parameters
    distances = zeros(1, nReal);
    velocities = zeros(1, nReal);
    for i = 1:nReal
        delay = real(lambdas(i))/(2*pi*SCS);
        fd = imag(lambdas(i))*SCS/(2*pi);
        distances(i) = -delay*c;
        velocities(i) = fd*c/fc;
    end
    results = [velocities', distances'];
end

function Q = getQTransformationMatrix(N)
    if mod(N, 2) == 0   % N even
        I = eye(N/2);
        J = flipud(I);
        Q = (1/sqrt(2))*[I, 1i*I; J, -1i*J];
    else                % N odd
        I = eye((N-1)/2);
        J = flipud(I);
        Z = zeros((N-1)/2, 1);
        Q = (1/sqrt(2))*[I, Z, 1i*I; Z', sqrt(2), Z'; J, Z, -1i*J];
    end
end

function ix = detectEigJumpBackwards(eigs, jumpThreshold, n)
    % obtain a threshold for eigenvalues based on samples that are from the
    % noise subspace and compare eigenvalues that could be from the signal
    % subspace to this threshold
    eigsNoise = mean(eigs(n+1:end));
    ix = 1;
    for i = 2:n
        if eigs(i) > eigsNoise*jumpThreshold
            ix = ix + 1;
        end
    end
end