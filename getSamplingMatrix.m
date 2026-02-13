function [PhiN, PhiM, nDiff] = getSamplingMatrix(Nmax, Mmax, Nsamples, Msamples)
    PhiN = zeros(Nsamples, Nmax);
    PhiM = zeros(Msamples, Mmax);
    nIxs = sort(randperm(Nmax, Nsamples));
    mIxs = sort(randperm(Mmax, Msamples));
    for i = 1:length(nIxs)
        PhiN(i, nIxs(i)) = 1;
    end
    for i = 1:length(mIxs)
        PhiM(i, mIxs(i)) = 1;
    end
    nDiff = nIxs(end) - nIxs(1);
end