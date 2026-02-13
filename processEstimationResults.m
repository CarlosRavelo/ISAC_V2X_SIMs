function [error, nDetected] = processEstimationResults(trueDistances, trueVelocities, results)
    nScatterers = length(trueDistances);
    nDetected = 0;

    scenarioTrue = [trueVelocities', trueDistances'];
    % separate target values
    targetTrue = scenarioTrue(end, :);
    % search estimation closest to target
    dToTarget = zeros(nDetected,1);
    for i = 1:size(results,1)
        dToTarget(i) = sqrt(sum((results(i, :)-targetTrue).^2));
    end
    [~, iClosest] = min(dToTarget);
    error = abs(results(iClosest, :)-targetTrue);
    [~, iToTarget] = sort(dToTarget);
    
    if iToTarget > 1
        for i = iToTarget(1:2)
            estTarget = results(i, :);
            dToEst = sqrt(sum((scenarioTrue-repmat(estTarget, [nScatterers, 1])).*(scenarioTrue-repmat(estTarget, [nScatterers, 1])), 2));
            [minDist, minToEst] = min(dToEst);
            if (minToEst == nScatterers) && (minDist < 2*norm(error))
                nDetected = 1;
                break
            end
        end
    else
        estTarget = results(iToTarget(1), :);
        dToEst = sqrt(sum((scenarioTrue-repmat(estTarget, [nScatterers, 1])).*(scenarioTrue-repmat(estTarget, [nScatterers, 1])), 2));
        [minDist, minToEst] = min(dToEst);
        if (minToEst == nScatterers) && (minDist < 2*norm(error))
                nDetected = 1;
        end
    end
end