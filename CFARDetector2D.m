function [peaks, th, cellsToTest] = CFARDetector2D(p, configStruct)
    % wrapper for the CFAR2D implemented in Matlab
    % p                     :       2D periodogram for peak detection (full or relevant section)
    % configStruct          :       structure with configuration parameters
    %   trainingBandSize    :       size of the training band
    %   guardBandSize       :       size of the guard band
    %   pfa                 :       probability of false alarm
    trainingBandSize = configStruct.trainingBand;
    guardBandSize = configStruct.guardBand;
    pfa = configStruct.pfa;
    
    detector = phased.CFARDetector2D('TrainingBandSize', trainingBandSize, ...
        'ThresholdFactor','Auto', 'GuardBandSize', guardBandSize, ...
        'ThresholdOutputPort', true, 'ProbabilityFalseAlarm', pfa);
    [N, M] = size(p);
    Ng = max(detector.GuardBandSize);
    Nt = max(detector.TrainingBandSize);
    cellsToTest = [];
    colstart = Nt + Ng + 1;
    rowstart = colstart;
    colend = M - (Nt + Ng);
    rowend = N - (Nt + Ng);

    for n = rowstart:rowend
        for m = colstart:colend
            cellsToTest = [cellsToTest, [n;m]];
        end
    end
    [peaks, th] = detector(p, cellsToTest);
end