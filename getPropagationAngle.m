function [az, el] = getPropagationAngle(posTx, posScatterer)
    % identify if the coordinates are in 2D or 3D
    if length(posTx) == 2
        el = 0;
        az = atand((posScatterer(1) - posTx(1))/(posScatterer(2) - posTx(2)));
    elseif length(posTx) == 3
        az = atand((posScatterer(1) - posTx(1))/(posScatterer(2) - posTx(2)));
        dz = posScatterer(3) - posTx(3);
        dxy = sqrt((posScatterer(1) - posTx(1))^2 + (posScatterer(2) - posTx(2))^2);
        el = atand(dz/dxy);
    else
        error('Wrong dimmensions for position vector')
    end
end