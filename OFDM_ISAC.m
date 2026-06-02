classdef OFDM_ISAC < handle
    properties
        grid
        fc
        SCS
    end

    methods
        function obj = OFDM_ISAC(grid, fc, SCS)
            %OFDM_ISAC creates an OFDM_ISAC object
            % obj = OFDM_ISAC(grid, fc, SCS) create an OFDM_ISAC object
            %
            %   Inputs:
            %       grid        :   symbol grid
            %       fc          :   carrier frequency
            %       SCS         :   subcarrier spacing
            obj.grid = grid;
            obj.fc = fc;
            obj.SCS = SCS;
        end

        function grid = getOFDMGrid(obj, trxObject, posScatterers, velScatterers, ...
                rcsScatterers)
            %GETOFDMGRID obtains the resulting OFDM grid 
            % GETOFDMGRID(trxObject, posScatterers, velScatterers,
            % rcsScatterers) estimates and aggregates the multipath of each
            % scatterer passed
            % 
            % Inputs:
            %   trxObjec        : object of the class ISAC_TRX
            %   posScatterers   : position of scatterers. Dimmensions are
            %                     number_of_scatterers x 2 with [x, y]
            %   velScatterers   : velocity of scatterers. Same dimmensions
            %                     with [vel_x, vel_y]
            %   rcsScatterers   : radar cross-section of scatterers in dBsm

            % get grid dimensions
            [N, M] = size(obj.grid);
            % initialize empty grid
            grid = zeros(N, M);
            for i = 1:length(rcsScatterers)
                
                % obtain path loss
                pathLossI = obj.getMonoRadarPathLoss(trxObject, ...
                                                    posScatterers(i, :), ...
                                                    rcsScatterers(i));
                [delay, doppler] = obj.getDelayDoppler(trxObject, ...
                    posScatterers(i, :), velScatterers(i, :));

                % get transmission angle
                [az, elev] = obj.getPropagationAngle(trxObject, posScatterers(i, :));
    
                % get the steering vector
                steeringVector = trxObject.steeringObject(obj.fc, [az; elev]);
    
                % get the Trx gain assumint Tx gain = Rx gain
                trxGain = trxObject.precoding'*(steeringVector*steeringVector')*trxObject.precoding;
    
                % get delay and doppler vectors
                delayVector = exp(-1j*2*pi*delay*(0:N-1)'*obj.SCS);
                dopplerVector = exp(1j*2*pi*(0:M-1)*doppler/obj.SCS);

                % obtain contribution
                pathI = sqrt(pathLossI*trxGain*trxObject.txPower/2)*delayVector*dopplerVector;    
                
                % add contribution
                grid = grid + pathI;
            end
        end

        function noisyGrid = getNoiseGrid(obj, effBW)
            % GETNOISEGRID Obtains the OFMD grid with the white gaussian
            % noise with the noise power determined by the effective
            % bandwidth
            % noisyGrid = GETNOISEGRID(effBW) returns the noise grid
            % Inputs:
            %       effBW       : effective bandwidth (Hz)
            % Output:
            %       noisyGrid   : ofdm grid with white gaussian noise
            T = 290;
            noisePower = physconst('Boltzmann')*T*effBW;
            [N, M] = size(obj.grid);
            noisyGrid = sqrt(noisePower/2)*complex(randn(N, M), randn(N, M));
        end
    end

    methods (Access=private)
        function pathLoss = getMonoRadarPathLoss(obj, trxObj, posScatterer, rcsScatterer)
            c0 = physconst('LightSpeed');
            rcsScatterer = 10^(rcsScatterer/10);
            d = sqrt(sum((trxObj.position - posScatterer).^2));
            pathLoss = c0^2*rcsScatterer/((4*pi)^3*d^4*obj.fc^2);
        end

        function [delay, doppler] = getDelayDoppler(obj, trxObject, ...
                                                    posScatterer, velScatterer)
            c0 = physconst('LightSpeed');

            % calculate the delay
            distance = norm(trxObject.position - posScatterer);
            delay = 2*distance/c0;
        
            % calculate the doppler
            velRad = obj.getRadialVelocity(trxObject, posScatterer, velScatterer);
            doppler = 2*velRad*obj.fc/c0;

            
        end
    end

    methods(Access=private, Static)
        function radialVelocity = getRadialVelocity(trxObject, ...
                                            posScatterer, velScatterer)
            % 1. obtain the position vector
            dVector = posScatterer - trxObject.position;
            dVector = dVector/norm(dVector);
            % 2. project source velocity into vector
            vRelTx = dot(dVector, trxObject.velocity);   
            % 3. project scatterer velocity into vector
            velRelScat = dot(dVector, velScatterer);
            % 4. obtain radial velocity
            radialVelocity = vRelTx - velRelScat;
        end

        function [az, elev] = getPropagationAngle(trxObj, posScatterer)
            % identify if the coordinates are in 2D or 3D
            posTx = trxObj.position;
            if length(posTx) == 2
                elev = 0;
                az = atand((posScatterer(1) - posTx(1))/(posScatterer(2) - posTx(2)));
            elseif length(posTx) == 3
                az = atand((posScatterer(1) - posTx(1))/(posScatterer(2) - posTx(2)));
                dz = posScatterer(3) - posTx(3);
                dxy = sqrt((posScatterer(1) - posTx(1))^2 + (posScatterer(2) - posTx(2))^2);
                elev = atand(dz/dxy);
            else
                error('Wrong dimmensions for position vector')
            end
        end
    end
end