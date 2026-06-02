% Example of usage of estimators for a simulation. It can be used as a
% guideline to adapt the code to different scenarios and models.
%% Parameter configuration

rng(1)

% Basic variables
c0  = physconst("LightSpeed");
T   = 290;                          % temperature in Kelvins (for noise power)

% Signal parameters
numerology  = 3;
slotSymbols = 14;                   % number of symbols per slot (14 for normal CP)
fc          = 28e9;                 % carrier frequency
lambda      = c0 / fc;              % wavelength
txPowerdBm  = 23;
txPower     = 10^(txPowerdBm/10)/1000;
BW          = 200;                  % total bandwidth in Hz
SCS         = 120e3;                % subcarrier spacing in Hz
PRBs        = 132;                  % total number of Physical Resource Blocks
Nmax        = PRBs*12;              % max number of subcarriers
numSlots    = 10*2^numerology;      % number of slots
Mmax        = slotSymbols*numSlots; % max number of OFDM symbols
emptyGrid   = zeros(Nmax, Mmax);    % empty OFDM grid
% ----------------
ofdmISAC    = OFDM_ISAC(emptyGrid, ...
                        fc, SCS);   % create the OFDM_ISAC object
% ----------------


% Estimation limits
dMax        = 70;                   % maximum estimation distance (m)
vMax        = 30;                   % maximum estimation radial velocity (m/s)

% ------ Estimation Method Parameter Configuration ------

% #### Periodogram ####
deltaN          = 12;                   % space between sampled subcarriers
deltaM          = 5;                    % space between sampled symbols
nSubcarriers    = 120;                  % used subcarriers
mSymbols        = 28;                   % used symbols
window          = 'hamming';            % window used for the periodogram
% ----------------
periodogram2D   = Periodogram2D(nSubcarriers, mSymbols, deltaN, ...
                                deltaM, window);    % Periodogram2D object
% ----------------
% cfar parameters
pfa             = 1e-6;                 % probability of false alarm
trainingBand    = [4, 4];               % size of training band
guardBand       = [1, 1];               % size of guard band
cfarDetector    = phased.CFARDetector2D('TrainingBandSize', trainingBand, ...
    'ThresholdFactor', 'Auto', 'GuardBandSize', guardBand, ...
    'ThresholdOutputPort', true, 'ProbabilityFalseAlarm', pfa);
% cluster detection parameters
epsilon         = 5;                    % cluster search radius
minpts          = 1;                    % minimum number of neighbors required
                                        % for core point

% #### MUSIC2D ####
nSubcarriers    = 120;                  % used subcarriers
mSymbols        = 28;                   % used symbols
lagN            = 8;                    % size of subarrays for rows
lagM            = 8;                    % size of subarrays for columns
dSearchSpace    = linspace(0, dMax, 100);       % distance search space
vSearchSpace    = linspace(-vMax, vMax, 100);   % velocity search space
FBLn            = 15;                   % rows of 2D full-back averaged matrix
FBLm            = 15;                   % columns of 2D full-back averaged matrix
musicOrderThrehsold = 100;                  % threshold for model order detection
decimationFactor= 3;                    % decimation factor for peak matching
peakProminence  = 1.5;                  % prominence for peak search
% ----------------
music2D         = MUSIC2D(nSubcarriers, mSymbols, lagN, lagM, ...
                    dSearchSpace, vSearchSpace, peakProminence, ...
                    FBLn, FBLm, musicOrderThrehsold, decimationFactor);
% ----------------

% #### ESPRIT2D ####
nSubcarriers    = 120;                  % used subcarriers
mSymbols        = 28;                   % used symbols
lagN            = 8;                    % rows in each snapshot
lagM            = 8;                    % columns in each snapshot
nIn             = 20;                   % rows in snapshot subsamples
mIn             = 8;                    % columns in snapshot subsamples
espritOrderThreshold = 5;               % threshold for model order detection
orderMax        = 13;                   % maximum model order
% ----------------
esprit2D = ESPRIT2D(nSubcarriers, mSymbols, lagN, lagM, nIn, mIn, ...
    espritOrderThreshold, orderMax);    % create the ESPRIT2D object

%% Scenario
% Ego vehicle parameters
egoVehiclePos   = [0, 0];
egoVehicleVel   = [0, 15];              % m/s
nElements       = 8;                    % number of antenna elements
trxAntenna      = phased.ULA(nElements, lambda/2);
trxObject       = ISAC_TRX(egoVehiclePos, egoVehicleVel, trxAntenna, ...
                    txPower);
trxObject.setPrecoding(fc, [0; 0]);

% Target parameters
targetRCSdB     = 10;                   % dBsm
targetRCS       = 10^(targetRCSdB/10);
targetPos       = [0, 30];              % m
targetVel       = [0, 0];               % m/s

% Street parameters
street.streetWidth          = 5;        % m
street.streetLength         = 70;       % m
street.sidewalkWidth        = 3;        % m

% Scatterer parameters
numScatterers               = 10;
rcsParams.rcsMean           = 5;        % dBsm
rcsParams.rcsVar            = 2;        % dBsm2

scatterers = placeSidewalkScatterers(numScatterers, street, rcsParams);

%% Perform estimatiom
method = {'periodogram', 'music2D', 'esprit2D'};

% contributions of scatterers
channelGridScatts = ofdmISAC.getOFDMGrid(trxObject, scatterers.positions, ...
    scatterers.velocities, scatterers.RCSs);
% contribution of target
channelGridTarget = ofdmISAC.getOFDMGrid(trxObject, targetPos, targetVel, ...
                    targetRCS);

channelGrid = channelGridScatts + channelGridTarget;

% select method to use
methodIx = 3;

switch method{methodIx}
    case 'periodogram'
        effBW           = periodogram2D.getEffectiveBandwidth(ofdmISAC);
        noisyGrid       = ofdmISAC.getNoiseGrid(effBW);
        ofdmISAC.grid   = channelGrid + noisyGrid;
        periodogram     = periodogram2D.getPeriodogram2D(ofdmISAC);
        results         = periodogram2D.estimateTargets(periodogram, ...
                            cfarDetector, epsilon, minpts);

    case 'music2D'
        effBW           = music2D.getEffectiveBandwidth(ofdmISAC);
        noisyGrid       = ofdmISAC.getNoiseGrid(effBW);
        ofdmISAC.grid   = channelGrid + noisyGrid;
        results         = music2D.getMUSIC2DEstimation(ofdmISAC);

    case 'esprit2D'
        effBW           = esprit2D.getEffectiveBandwidth(ofdmISAC);
        noisyGrid       = ofdmISAC.getNoiseGrid(effBW);
        ofdmISAC.grid   = channelGrid + noisyGrid;
        results         = esprit2D.getESPRIT2DEstimation(ofdmISAC);
    otherwise
        error('Invalid method selected')
end




