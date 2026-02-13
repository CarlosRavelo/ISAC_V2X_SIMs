function [Z, indexes] = fomp2D(Y, A, A_t, l)
    R = Y;
    [~, N] = size(A);
    [M, ~] = size(A_t);
    F = ones(N, M);
    indexes = zeros(l, 2);
    [Ns, Ms] = size(Y);
    Bmatrix = zeros(Ns, Ms, l);
    for i = 1:l
        % 1. obtain the biggest contribution
        Zi = A'*R*A_t'.*F;
        [x, n] = max(Zi); [~, m] = max(x);
        nI = n(m);
        mI = m;
        % remove contribution for future iterations
        F(nI, mI) = 0;
        % add indexes to record
        indexes(i, :) = [nI, mI];
        
        % 2. get the weights for the selected atoms
        % obtain H and f simultaneously
        H = zeros(i);
        f = zeros(i, 1);
        Bmatrix(:, :, i) = A(:, indexes(i, 1))*A_t(indexes(i, 2), :);
        for j = 1:i          
            f(j) = inner(Y, Bmatrix(:, :, j));
            for k = 1:i
                H(j, k) = inner(Bmatrix(:, :, j), Bmatrix(:, :, k));
            end
        end
        uWeights = H\f;

        % 3. update residue
        R = Y;
        for j = 1:i
            Bj = A(:, indexes(j, 1))*A_t(indexes(j, 2), :);
            R = R - uWeights(j)*Bj;
        end

        % 4. check residue error
        norm(Y-R, 'fro');
    end
    % move the final weights to the adequate positions in the estimation
    Z = zeros(N, M);
    for j = 1:i
        Z(indexes(j, 1), indexes(j, 2)) = uWeights(j);
    end
end

function result = inner(A, B)
    result = sum(A(:).*B(:));
end