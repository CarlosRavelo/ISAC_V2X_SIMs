function uniquePeaks = getClusterPeaks(p, peaks, testedCells, epsilon, minpts)
    % detect clusters of peaks from CFAR and take only the dominant point
    % for each cluster
    % p                 :       data over which peaks were detected
    % peaks             :       cells where peaks were detected
    % testedCells       :       cells where the peak search worked
    % epsilon           :       maximum distance between clusters
    % minpts            :       minimum points for a cluster
    peakCells = testedCells(:, peaks)';
    clusters = dbscan(peakCells, epsilon, minpts);
    uniqueClusters = unique(clusters);
    peakValues = p(sub2ind(size(p), peakCells(:, 1), peakCells(:, 2)));
    uniquePeaks = zeros(length(uniqueClusters), 2);

    for i = 1:length(uniqueClusters)
        clusterIndices = find(clusters == uniqueClusters(i));
        clusterPeakValues = peakValues(clusterIndices);
        [~, idx] = max(clusterPeakValues);
        uniquePeaks(i, :) = peakCells(clusterIndices(idx), :);
    end
end