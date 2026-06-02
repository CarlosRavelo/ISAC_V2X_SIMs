classdef ESPRIT2D
    properties
        nSubcarriers(1, 1)          {mustBeInteger, mustBeNonnegative}
        mSymbols(1, 1)              {mustBeInteger, mustBeNonnegative}
        lagN(1, 1)                  {mustBeInteger, mustBeNonnegative}
        lagM(1, 1)                  {mustBeInteger, mustBeNonnegative}
        nIn(1, 1)                   {mustBeInteger, mustBeNonnegative}
        mIn(1, 1)                   {mustBeInteger, mustBeNonnegative}
        modelOrderThreshold(1, 1)   {mustBeNonnegative}
        orderMax(1, 1)              {mustBeInteger, mustBeNonnegative}
    end

    methods
        function obj = ESPRIT2D(nSubcarriers, mSymbols, lagN, lagM, nIn, mIn, ...
                modelOrderThreshold, orderMax)
            % ESPRIT2D Creates a ESPRIT2D object
            % obj = ESPRIT2D(...) creates the ESPRIT2D object with the
            % specified parameters
            %
            % Inputs:
            %       nSubcarriers        : number of subcarriers to be used
            %       mSymblols           : number of symbols to be used
            %       lagN                : rows in each snapshot
            %       lagM                : columns in each snapshot
            %       nIn                 : rows in snapshot subsamples
            %       mIn                 : columns in snapshot subsamples
            %       modelOrderThreshold : threshold for model order
            %                             estimation
            %       orderMax            : upper threshold for model order
            % Output:
            %       obj: ESPRIT2D object
        
            obj.nSubcarriers = nSubcarriers;
            obj.mSymbols = mSymbols;
            obj.lagN = lagN;
            obj.lagM = lagM;
            obj.nIn = nIn;
            obj.mIn = mIn;
            obj.modelOrderThreshold = modelOrderThreshold;
            obj.orderMax = orderMax;
        end

        function results = getESPRIT2DEstimation(obj, ofdm)
            % GETESPRIT2DESTIMATION Estimates distance and velocity through
            % 2D ESPRIT
            % results = GETESPRIT2DESTIMATION(ofdm) returns the
            % distance and velocity of scatterers in scenario. 
            % 
            % Inputs:
            %       ofdm                : received grid on sensing receiver
            % Output:
            %       results             : Ns x 2 matrix with distance and
            %                             velocity estimation for Ns
            %                             scatterers detected. results(i,
            %                             :) = [distance(Ni), velocity(Ni)]
            arguments
                obj
                ofdm OFDM_ISAC
            end

            c0 = physconst('LightSpeed');
            
            rx = ofdm.grid(1:obj.nSubcarriers, 1:obj.mSymbols);
            [N, M] = size(rx);
            kRows = N - obj.lagN + 1;
            kColumns = M - obj.lagM + 1;

            % get auxiliary matrices
            % unitary transformation matrices
            Qn = obj.getQTransformationMatrix(obj.nIn);
            Qn1 = obj.getQTransformationMatrix(obj.nIn - 1);
            Qm = obj.getQTransformationMatrix(obj.mIn);
            Qm1 = obj.getQTransformationMatrix(obj.mIn - 1);

            % selection matrices
            Js1N = [eye(obj.nIn-1), zeros(obj.nIn-1, 1)];
            Js2N = [zeros(obj.nIn-1, 1), eye(obj.nIn-1)];
            Js2M = [zeros(obj.mIn-1, 1), eye(obj.mIn-1)];

            K1 = real(Qn1'*Js1N*Qn);
            K2 = imag(Qn1'*Js2N*Qn);
            K3 = real(Qm1'*Js2M*Qm);
            K4 = imag(Qm1'*Js2M*Qm);

            Kx1 = kron(eye(obj.mIn), K1);
            Kx2 = kron(eye(obj.mIn), K2);
            Ky1 = kron(K3, eye(obj.nIn));
            Ky2 = kron(K4, eye(obj.nIn));

            % snapshots
            Y = zeros(obj.nIn*obj.mIn, kRows*kColumns);
            index = 1;

            nSubsamples = obj.lagN - obj.nIn + 1;
            mSubsamples = obj.lagM - obj.mIn + 1;

            for i = 1:kRows
                for k = 1:kColumns
                    rxI = rx(i:i+obj.lagN-1, k:k+obj.lagM-1);
                    for ii = 1:nSubsamples
                        for kk = 1:mSubsamples
                            rxIn = rxI(ii:ii+obj.nIn-1, kk:kk+obj.mIn-1);
                            Yii = reshape(Qn'*rxIn*conj(Qm), obj.nIn*obj.mIn, 1);
                            Y(:, index) = Y(:, index) + Yii;
                        end
                    end
                end
            end

            Ymat = [real(Y), imag(Y)];

            % obtain the SVD
            [U, P, ~] = svd(Ymat);
            % estimate model order
            order = obj.estimateModelOrder(diag(P), obj.modelOrderThreshold, ...
                obj.orderMax);
            % signal subspace
            Us = U(:, 1:order);

            Psix = (Kx1*Us)\(Kx2*Us);
            Psiy = (Ky1*Us)\(Ky2*Us);
            PsiMatrix = Psix + 1j*Psiy;
            lambdas = eig(PsiMatrix);

            distances = -real(lambdas(1:order))*c0/(2*pi*ofdm.SCS);
            velocities = imag(lambdas(1:order))*ofdm.SCS*c0/(2*pi*ofdm.fc);

            results = [distances, velocities];
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
            effBW = obj.nSubcarriers * ofdmObj.SCS;
        end
    end

    methods (Access=private, Static)
        function Q = getQTransformationMatrix(N)
            % GETQTRANSFORMATIONMATRIX Forms unitary transformation matrix
            % Q = GETQTRANSFORMATIONMATRIX(N) returns the unitary
            % transformation matrix Q
            % Input:
            %       N       : target dimension for Q
            % Output:
            %       Q       : unitary trasnformation matrix of size NxN
            if mod(N, 2) == 0   % even
                I = eye(N/2);
                J = flipud(I);
                Q = (1/sqrt(2))*[I, 1j*I; J, -1j*J];
            else                % odd
                I = eye((N-1)/2);
                J = flipud(I);
                z = zeros((N-1)/2, 1);
                Q = (1/sqrt(2))*[I, z, 1j*I; z', sqrt(2), z'; J, z, -1j*J];
            end
        end

        function ix = estimateModelOrder(eigs, threshold, maxOrder)
            % ESTIMATEMODELORDER Estimates model order based on the
            % eigenvalues of the SVD of the correlation matrix
            % ix = ESTIMATEMODELORDER(...) returns model order estimation
            % Input:
            %       eigs        : eigenvalues of SVD of correlation matrix
            %       threshold   : threshold for order estimation
            %       maxOrder    : upper limit for model order
            % Output:
            %       ix          : model order estimation
            eigsNoise = mean(eigs(maxOrder+1:end));
            ix = 1;
            for i = 1:maxOrder
                if eigs(i) > eigsNoise*threshold
                    ix = ix + 1;
                else
                    break
                end
            end
        end
    end
end