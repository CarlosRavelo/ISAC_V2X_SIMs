function noisePower = estimateNoisePower(T, BW)
    % estimate the noise power for a given bandwidth and ambient
    % temperature
    noisePower = physconst('Boltzmann')*T*BW;
end