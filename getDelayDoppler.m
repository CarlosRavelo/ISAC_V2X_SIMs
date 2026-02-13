function [delay, doppler, distance, velRad] = getDelayDoppler(posTx, posScatterer, velTx, velScatterer, fc)
    % obtain the delay and doppler, given the positions of the source and
    % scatterer on a 2D Plane. 
    C = physconst('LightSpeed');

    % firstly obtain the delay
    distance = norm(posTx - posScatterer);
    delay = 2*distance/C;

    % calculate the doppler
    velRad = getRadialVelocity(posTx, posScatterer, velTx, velScatterer);
    doppler = 2*velRad*fc/C;
end