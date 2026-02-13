%% Process batch of simulations for testing the parameter estimation methods
rng(1)

% Basic variables
C = physconst('LightSpeed');
T = 290;

% Configure the signal parameters
fc = 24e9;
lambda = C/fc;
txPowerdBm = 23;
txPower = 10^(txPowerdBm/10)/1000;
BW = 200e6;
SCS = 120e3;
PRBs = 132;
Nmax = PRBs*12;
numSlots = 10*8;
Mmax = 14*numSlots;
dMax = 70;
vMax = 30;

% Configure the processing method parameters
% periodogram
% first define the subsampling
deltaN = 12;
deltaM = 5;
Nuse = 64;
Muse = 32;
nIxUse = 1:deltaN:(Nuse)*deltaN;
mIxUse = 1:deltaM:(Muse)*deltaM;
periodogramParams = struct;
periodogramParams.N = Nmax;
periodogramParams.M = Mmax;
periodogramParams.EffBW = SCS*periodogramParams.N;
periodogramParams.Window = 'hamming';
periodogramParams.dMax = dMax;
periodogramParams.vMax = vMax;
periodogramParams.deltaN = deltaN;
periodogramParams.deltaM = deltaM;
cfarParams = struct;
cfarParams.pfa = 1e-6;
cfarParams.trainingBand = [4, 4];               % for reduced resources
% cfarParams.trainingBand = [8, 8];
cfarParams.guardBand = [1, 1];


% music
musicParams = struct;
musicParams.N = 64;
musicParams.M = 32;
musicParams.corrMatMethod = 1;
musicParams.Ln = 15;
musicParams.Lm = 15;
musicParams.dSearchSpace = linspace(0, 75, 100);
musicParams.vSearchSpace = linspace(-30, 30, 100);
peakSearchParams = struct;
peakSearchParams.peakProminence = 1.5;
peakSearchParams.FBLn = 15;
peakSearchParams.FBLm = 15;

% esprit
espritParams = struct;
espritParams.N = 64;
espritParams.M = 32;
espritParams.Nlag = 24;
espritParams.Mlag = 12;
espritParams.Nin = 20;
espritParams.Min = 8;

% omp
ompParams = struct;
ompParams.N = 256;
ompParams.M = 128;


% Configure the batch processing parameters
baseFilename = 'scattererParams';
specificFilename = [3, 6, 9, 12];
allMethods = {'periodogram', 'music', 'esprit', 'omp'};


% Ego vehicle parameters
egoVehiclePos = [0, 0];
egoVehicleVel = [0, 15];
nAntennas = 8;
trxAntenna = phased.ULA(nAntennas, lambda/2);
trxSteerer = phased.SteeringVector('SensorArray', trxAntenna);
trxPrecodingVector = trxSteerer(fc, [0; 0]);

% Configure the target parameters
targetRCSdB = 10;
targetDistanceRange = 10:5:70;
targetVel = [0, 10];

% Process each scenario
for ix = 1:4
    selectedMethod = allMethods{ix};
    for i = specificFilename
        % 1. Import the scenario locations
        fileI = load(strcat(baseFilename, int2str(i), '.mat'));
        scenarioI = fileI.scattererParams;
        % errorsI = zeros(length(scenarioI), length(targetDistanceRange), 3);
        errorsI = zeros(1, length(targetDistanceRange), 3);
        for j = 1:length(scenarioI)
        % for j = 1:1
            % obtain the OFDM grid from all scatterers
            ofdmGrid = zeros(Nmax, Mmax);
            scatterersJ = scenarioI{j};
            scatterersDistances = zeros(1, length(scatterersJ.RCSs));
            scatterersRadVelocity = zeros(1, length(scatterersJ.RCSs));
            for jj = 1:length(scatterersJ.Positions)
                [scattererChannelJJ, pathLossJJ, scatterersDistances(jj), scatterersRadVelocity(jj)] = getChannelOFDMGrid(egoVehiclePos, ...
                                                                scatterersJ.Positions(:, jj)', egoVehicleVel, [0, 0], ...
                                                                scatterersJ.RCSs(jj), fc, trxAntenna, ...
                                                                trxPrecodingVector, Nmax, Mmax, SCS);
                ofdmGrid = ofdmGrid + scattererChannelJJ;
            end
            % obtain and process the target contribution for every distance
            errorsPerDistance = zeros(2, i, length(targetDistanceRange));
            detectionsPerDistance = zeros(1, length(targetDistanceRange));
            for d = 1:length(targetDistanceRange)
                [targetChannel, targetPathLoss, ~, targetVelocity] = getChannelOFDMGrid(egoVehiclePos, [0, targetDistanceRange(d)], ...
                                        egoVehicleVel, targetVel, ...
                                        targetRCSdB, fc, trxAntenna, ...
                                        trxPrecodingVector, Nmax, Mmax, SCS);
                ofdmGridAll = sqrt(txPower/2)*(ofdmGrid+targetChannel);
                % process the scenario
                switch selectedMethod
                    case 'periodogram'
                        % estimate the noise and add it to the signal
                        nPower = estimateNoisePower(T, periodogramParams.EffBW);
                        ofdmGridAll = ofdmGridAll + sqrt(nPower/2)*complex(randn(Nmax, Mmax), randn(Nmax, Mmax));
    
                        % subsample the periodogram
                        ofdmGridAll = ofdmGridAll(nIxUse, mIxUse);
    
                        % call the periodogram method
                        [periodogram, dAxis, vAxis] = getPeriodogram(ofdmGridAll, fc, SCS, periodogramParams);
                        periodogramDB = 10*log10(periodogram);
                        % plot if I am debugging
                        % imagesc(vAxis, dAxis, periodogramDB)
    
                        % run the CFAR detection method
                        [peaks, th, testedCells] = CFARDetector2D(periodogram, cfarParams);
                        
                        % detect clusters of peaks
                        epsilon = 5;
                        minpts = 1;
                        uniquePeaks = getClusterPeaks(periodogramDB, peaks, testedCells, epsilon, minpts);
                        results = [vAxis(uniquePeaks(:, 2))', dAxis(uniquePeaks(:, 1))'];
                    case 'music'
                        % add the noise to the signal
                        nPower = estimateNoisePower(T, musicParams.N*SCS);
                        noiseGrid = sqrt(nPower/2)*complex(randn(musicParams.N, musicParams.M), ...
                                                           randn(musicParams.N, musicParams.M));
                        ofdmGridAll = ofdmGridAll(1:musicParams.N, 1:musicParams.M) + noiseGrid;
    
                        % call the 1D method on both dimmensions
                        musicParams.Order = i+1;
                        [psn, psm, nEst] = music1D(ofdmGridAll, SCS, fc, musicParams);
                        % transform to dB for better peak estimation
                        Psn = pow2db(psn/max(psn));
                        Psm = pow2db(psm/max(psm));
                        % figure
                        % plot(musicParams.dSearchSpace, Psn);
                        % title('Pseudospectrum distance')
                        % figure 
                        % plot(musicParams.vSearchSpace, Psm);
                        % title('Pseudospectrum velocity')
    
                        % search for the peaks
                        peakSearchParams.nOrder = nEst;
                        peakSearchParams.dSpace = musicParams.dSearchSpace;
                        peakSearchParams.vSpace = musicParams.vSearchSpace;
                        
                        results = getMUSICPeaks(ofdmGridAll, Psn, Psm, SCS, fc, peakSearchParams);
                        
                    case 'esprit'
                        nPower = estimateNoisePower(T, espritParams.N*SCS);
                        noiseGrid = sqrt(nPower/2)*complex(randn(espritParams.N, espritParams.M), ...
                                                            randn(espritParams.N, espritParams.M));
                        rx = ofdmGridAll(1:espritParams.N, 1:espritParams.M) + noiseGrid;
                        espritParams.Order = i + 1;
                        results = unitaryESPRIT2DFB(rx, SCS, fc, espritParams);
                    case 'omp'
                        % get the transformation matrix
                        Fn = dftmtx(Nmax);
                        Fm = dftmtx(Mmax);
    
                        % get the sampling matrix
                        [PhiN, PhiM, nDiff] = getSamplingMatrix(Nmax, Mmax, ompParams.N, ompParams.M);
    
                        % get the atoms
                        A = PhiN*Fn';
                        A_t = Fm'*PhiM';
    
                        % add noise depending on bandwidth sampled
                        nPower = estimateNoisePower(T, nDiff*SCS);
                        noiseGrid = sqrt(nPower/2)*complex(randn(Nmax, Mmax), randn(Nmax, Mmax));
                        ofdmGridAll = ofdmGridAll + noiseGrid;
    
                        % get the sparsely sampled signal
                        ySparse = PhiN*ofdmGridAll*PhiM';
                        
                        % correct this
                        ompParams.Order = i+1;
                        [Z, peaks] = fomp2D(ySparse, A, A_t, ompParams.Order);
                        Z = flipud(fftshift(Z, 2));
                       
    
                        % create axis and locate the peaks
                        dAxis = (0:Nmax-1)*C/(2*SCS*Nmax);
                        vAxis = (-Mmax/2:Mmax/2-1)*C*SCS/(2*fc*Mmax);
                        
                        peaksLoc = zeros(size(peaks));
                        vAxisHalf = length(vAxis)/2;
                        for ii = 1:size(peaks, 1)
                            peaksLoc(ii, 1) = length(dAxis) - peaks(ii, 1);
                            if peaks(ii, 2) < vAxisHalf
                                peaksLoc(ii, 2) = peaks(ii, 2) + vAxisHalf;
                            else
                                peaksLoc(ii, 2) = peaks(ii, 2) - vAxisHalf;
                            end
                        end
                        peaksLoc(peaksLoc==0) = 1;
    
                        results = [vAxis(peaksLoc(:, 2))', dAxis(peaksLoc(:, 1))'];
                end
                trueDistances = [scatterersDistances, targetDistanceRange(d)];
                % trueDistances = targetDistanceRange(d);
                trueVelocities = [scatterersRadVelocity, targetVelocity];
                % trueVelocities = targetVelocity;
                [errors, detectedTarget] = processEstimationResults(trueDistances, trueVelocities, results);
                errorsI(j, d, 1) = errors(1);
                errorsI(j, d, 2) = errors(2);
                errorsI(j, d, 3) = detectedTarget;
            end
        end
        % save the result to a file
        fileNameI = strcat('results', selectedMethod, int2str(i), '.mat');
        save(fileNameI, 'errorsI');
    end
end