classdef Periodogram2D
    properties
        nSubcarriers(1, 1)  {mustBeInteger, mustBeNonnegative}
        mSymbols(1, 1)  {mustBeInteger, mustBeNonnegative}
        deltaN(1, 1)  {mustBeInteger, mustBeNonnegative}
        deltaM(1, 1)  {mustBeInteger, mustBeNonnegative}
        window string 
    end

    methods
        function obj = Periodogram2D(nSubcarriers, mSymbols, ...
                deltaN, deltaM, window)
            %PERIODOGRAM2D Creates a Periodogram2D object
            %  obj = PERIODOGRAM2D(...) creates a Periodogram2D object with
            %  the specified parameters
            %
            %  Inputs:
            %       nSubcarriers: number of subcarriers to be used
            %       mSymblols   : number of symbols to be used
            %       deltaN      : spacing between used subcarriers
            %       deltaM      : spacing between used symbols
            %       window      : window to be used with the periodogram
            %                     must be one of hamming, blackman,
            %                     chebyshev or rect
            %   
            %
            %  Outputs:
            %       obj         : periodogram2D object
            obj.nSubcarriers = nSubcarriers;
            obj.mSymbols = mSymbols;
            obj.deltaN = deltaN;
            obj.deltaM = deltaM;
            obj.window = window;
            
        end
        function periodogram = getPeriodogram2D(obj, ofdm)
            %GETPERIODOGRAM2D calculates the 2D periodogram of the signal
            %provided in ofdm
            %  periodogram = GETPERIODOGRAM2D(...) calculates a 2D
            %  periodogram
            %  
            %  Inputs:
            %       obj         : Periodogram2D object
            %       ofdm        : object of the class OFDM_ISAC
            %
            %  Outputs:
            %       periodogram : periodogram structure including the 
            %                     distance and velocityaxis
            arguments
                obj
                ofdm OFDM_ISAC
            end

            grid = ofdm.grid;
            % use subsampling if indicated
            nUse = 1:obj.deltaN:obj.nSubcarriers*obj.deltaN;
            mUse = 1:obj.deltaM:obj.mSymbols*obj.deltaM;
            rx = grid(nUse, mUse);
            [N, M] = size(rx); 
            W = obj.getWindowMatrix(N, M);
            % get the 2D periodogram
            P = ifft(rx.*W, 2^nextpow2(N), 1);
            P = fft(P, 2^nextpow2(M), 2);
            P = abs(fftshift(P, 2).^2);
            periodogram.P = P;

            % obtain the periodogram's axis
            SCS = ofdm.SCS;
            fc = ofdm.fc;
            c0 = physconst('LightSpeed');
            [nP, mP] = size(P);
            distanceAxis = (0:nP-1)*c0/(2*SCS*nP)/obj.deltaN;
            velocityAxis = (-mP/2:mP/2-1)*c0*SCS/(2*fc*mP*obj.deltaM);
            periodogram.distanceAxis = distanceAxis;
            periodogram.velocityAxis = velocityAxis;
        end

        function targets = estimateTargets(obj, periodogram, detector, ...
                epsilon, minpts)
            %ESTIMATETARGETS detects the targets and estimates their
            %distance and radial velocity
            %  targets = ESTIMATETARGETS(...) detect peaks and estimate
            %  parameters
            %  
            %  Inputs:
            %       obj         : Periodogram2D object
            %       periodogram : periodogram structure
            %       detector    : phased.CFARDetector2D object
            %       epsilon     : search radius around a point
            %       minpts      : minimum number of neighbors required for
            %                     core point
            %
            %  Outputs:
            %       targets     : matrix of dimmensions NumberOfTargets x 2
            %                     where the first column is the distance 
            %                     and the second the velocity
            arguments
                obj
                periodogram struct
                detector phased.CFARDetector2D
                epsilon {mustBeInteger, mustBePositive}
                minpts {mustBeInteger, mustBePositive}
            end
            % apply CFAR to detect the reflectors
            detections = obj.getDetections(periodogram, detector);
            % get the indices of detected peaks
            peakCells = detections.cellsToTest(:, detections.peaks);
            peakCells = peakCells';
            % detect clusters within the identified peaks
            clusters = dbscan(peakCells, epsilon, minpts);
            uniqueClusters = unique(clusters);
            P = periodogram.P;
            % periodogram values at peak locations
            peakValues = P(sub2ind(size(P), peakCells(:, 1), peakCells(:, 2)));
            uniquePeaks = zeros(length(uniqueClusters), 2);
            
            % keep cluster index with higher value (cluster max peak)
            for i = 1:length(uniqueClusters)
                clusterIndices = find(clusters == uniqueClusters(i));
                clusterPeakValues = peakValues(clusterIndices);
                [~, idx] = max(clusterPeakValues);
                uniquePeaks(i, :) = peakCells(clusterIndices(idx), :);
            end
            % translate peaks to distance and velocity estimations
            targets = [periodogram.distanceAxis(uniquePeaks(:, 1))', ... 
                periodogram.velocityAxis(uniquePeaks(:, 2))'];
        end

        function effBW = getEffectiveBandwidth(obj, ofdmObj)
            % GETEFFECTIVEBANDWIDTH Gets the bandwidth occupied by the
            % sensing resources
            % effBW = GETEFFECTIVEBANDWIDTH(obj, ofdmObj) calculates total
            % bandwidth
            % Input:
            %       ofdmObj : object of the OFDM_ISAC class

            arguments
                obj 
                ofdmObj OFDM_ISAC 
            end
            effBW = obj.nSubcarriers * ofdmObj.SCS * obj.deltaN;
        end
    end

    methods (Access = private)
        function W = getWindowMatrix(obj, N, M)
            %GETWINDOWMATRIX obtains the 2D window matrix of dimmensions
            %NxM
            %  W = GETWINDOWMATRIX(...) gets the windowing matrix
            %  
            %  Inputs:
            %       obj         : Periodogram2D object
            %       N           : rows of the window
            %       M           : columns of the window
            %
            %  Outputs:
            %       W           : windowing matrix

            switch obj.window
                case 'hamming'
                    wN = hamming(N);
                    wM = hamming(M)';
                case 'blackman'
                    wN = blackman(N);
                    wM = blackman(M)';
                case 'chebyshev'
                    wN = chebwin(N);
                    wM = chebwin(M)';
                case 'rect'
                    wN = ones(N, 1);
                    wM = ones(1, M);                
                otherwise
                    error('Wrong window type');
            end
            W = wN * wM;
        end

        function detections = getDetections(obj, periodogram, detector)
            %GETDETECTIONS detect scatterers in 2D periodogram using CFAR
            %  detections = GETDETECTIONS(...) detects scatterers
            %  
            %  Inputs:
            %       obj         : Periodogram2D object
            %       periodogram : periodogram structure
            %       detector    : CFARDetector2D object
            %
            %  Outputs:
            %       detections  : structure containing the detected peaks,
            %                     tested cells and thresholds
            arguments
                obj
                periodogram struct
                detector    
            end
            [N, M] = size(periodogram.P);
            
            Ng = max(detector.GuardBandSize);
            Nt = max(detector.TrainingBandSize);

            startIdx = Nt + Ng + 1;
            rowStart = startIdx;
            colStart = startIdx;
            rowEnd = N - (Nt + Ng);
            colEnd = M - (Nt + Ng);

            % get cells to test
            [rowIdx, colIdx] = ndgrid(rowStart: rowEnd, colStart: colEnd);
            cellsToTest = [rowIdx(:)'; colIdx(:)'];

            % apply cfar detector
            [peaks, th] = detector(periodogram.P, cellsToTest);

            % construct output structure
            detections.peaks = peaks;
            detections.thresholds = th;
            detections.cellsToTest = cellsToTest;
        end
    end
end