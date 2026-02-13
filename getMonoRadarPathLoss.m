function pathLoss = getMonoRadarPathLoss(posTx, posScatterer, fc, rcs)
    C = physconst('LightSpeed');
    d = sqrt(sum((posScatterer - posTx).^2));
    pathLoss = C^2*rcs/((4*pi)^3*d^4*fc^2);
end