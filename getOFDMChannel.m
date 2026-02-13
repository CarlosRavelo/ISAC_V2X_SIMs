%% Obtain the OFDM channel given a set of scatterers. Should apply to scatterers and target
function ofdmChannel = getOFDMChannel(scattererParams, scenarioParams, transmitterParams)
    % get the OFDM channel grid for a set of scatterers with a given set of
    % scenario and transmitter parameters
    % Input Parameters:
    % scattererParams       :       structure with the parameters of the
    %                               scattererers
    % scenarioParams        :       structure with the scenario parameters
    % transmitterParams     :       structure with the transmitter
    %                               parameters

    % initialize the scatterer parameters
    scattererPos = scattererParams.Positions;
    scattererRCS = scattererParams.RCSs;
    scattererVelocities = scattererParams.Velocities;
    numScatterers = size(scattererPos, 2);

    % initialize the scenario parameters
    fc = scenarioParams.Fc;
    posTx = scenarioParams.TxPos;
    vTx = scenarioParams.TxVelocity;

    % initialize the transmitter parameters
    nSubcarriers = transmitterParams.N;
    mSymbols = transmitterParams.M;
    SCS = transmitterParams.SCS;
    trxAnt = transmitterParams.Antenna;
    trxPrecoding = transmitterParams.Precoding;
    txPow = transmitterParams.TxPower;

    % initialize OFDM Grid
    ofdmChannel = zeros(nSubcarriers, mSymbols);
    % get the channel contribution of each scatterer
    for i = 1:numScatterers
        % obtain path loss
        pathLoss = getMonoRadarPathLoss(posTx, scattererPos(i), ...
                                        fc, scattererRCS(i));
        % obtain delay and doppler
        [delay, doppler] = getDelayDoppler(posTx, scattererPos(i), vTx, ...
                                            scattererVelocities(i), fc);
        
        % obtain the trx gain
        trxGain = getTrxGain(posTx, scattererPos(i), trxPrecoding, trxAnt, fc);

        % obtain the total path loss including tx power, trx gain and
        % propagation losses
        effectivePathLoss = pathLoss*(txPow/trxAnt.numElements)*trxGain;

        % obtain the delay and doppler vectors
        delayVector = exp(-1j*2*pi*delay*(0:nSubcarriers-1)'*SCS);
        dopplerVector = exp(1j*2*pi*(0:mSymbols-1)*doppler/SCS);

        % add contribution to channel
        ofdmChannel = ofdmChannel + sqrt(effectivePathLoss)*delayVector*dopplerVector;
    end
    
end

function pathLoss = getMonoRadarPathLoss(posTx, posScatterer, fc, rcs)
    C = physconst('LightSpeed');
    d = sqrt(sum((posScatterer - posTx).^2));
    pathLoss = C^2*rcs/((4*pi)^3*d^4*fc^2);
end

function [delay, doppler] = getDelayDoppler(posTx, posScatterer, velTx, velScatterer, fc)
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

function velRad = getRadialVelocity(posTx, posScatterer, velTx, velScatterer)
    % 1. obtain the joining vector
    joiningVector = posScatterer - posTx;
    joiningVector = joiningVector/norm(joiningVector);
    % 2. project source velocity into joining vector
    vRelTx = dot(joiningVector, velTx);   
    % 3. project scatterer velocity into joining vector
    velRelScat = dot(joiningVector, velScatterer);
    % 4. obtain radial velocity
    velRad = vRelTx - velRelScat;
end

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

function trxGain = getTrxGain(posTx, posScatterer, trxPrecoding, trxAntenna, fc)
    % get the trx gain for a given precoding and antenna
    
    % first obtain the azimuth and elevation angle
    [az, el] = getPropagationAngle(posTx, posScatterer);
    % create the steering object and obtain the steering vector
    steeringObject = phased.SteeringVector('SensorArray', trxAntenna);
    steeringVector = steeringObject(fc, [az; el]);
    
    % obtain the trx gain
    trxGain = trxPrecoding'*(steeringVector*steeringVector')*trxPrecoding;
end

