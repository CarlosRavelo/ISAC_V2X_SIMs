function [ofdmChannel, pathLossPlusGain, distance, velRal] = getChannelOFDMGrid(posTx, posScatterer, velTx, velScatterer, ...
                                          rcsScatterer, fc, trxAntenna, ...
                                          precoding, N, M, SCS)
    % obtain path loss
    pathLoss = getMonoRadarPathLoss(posTx, posScatterer, fc, rcsScatterer);
    % obtain delay and doppler
    [delay, doppler, distance, velRal] = getDelayDoppler(posTx, posScatterer, velTx, velScatterer, fc);
    % get the angle of transmission and reception to get steerers
    [az, el] = getPropagationAngle(posTx, posScatterer);
    % get the steering vector
    steeringObject = phased.SteeringVector('SensorArray', trxAntenna);
    steeringVector = steeringObject(fc, [az; el]);
    % add tx and rx gain to the pathloss
    TrxGain = precoding'*(steeringVector*steeringVector')*precoding;
    pathLossPlusGain = pathLoss*TrxGain;
    % get a random phase shift from reflection
    % randomShift = rand()*2*pi;
    % create delay and doppler vectors
    delayVector = exp(-1j*2*pi*delay*(0:N-1)'*SCS);
    dopplerVector = exp(1j*2*pi*(0:M-1)*doppler/SCS);
    % obtain the channel in the OFDM grid
    %ofdmChannel = sqrt(pathLossPlusGain/2)*delayVector*dopplerVector.*exp(1j*randomShift);
    ofdmChannel = sqrt(pathLossPlusGain/2)*delayVector*dopplerVector;
end