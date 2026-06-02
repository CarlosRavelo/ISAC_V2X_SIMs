classdef ISAC_TRX < handle
    properties
        position
        velocity
        trxAntenna
        steeringObject
        precoding
        txPower
    end
    
    methods
        function obj = ISAC_TRX(position, velocity, trxAntenna, txPower)
            %ISAC_TRX Creates an ISAC_TRX object
            % obj = ISAC_TRX(position, velocity, trxAntenna, txPower
            % 
            % Inputs:
            %       position        : position of the Trx node
            %       velocity        : velocity of the Trx node
            %       trxAntenna      : ULA or URA object
            %       txPower         : transmit power in Watts
            obj.position = position;
            obj.velocity = velocity;
            obj.trxAntenna = trxAntenna;
            obj.steeringObject = phased.SteeringVector('SensorArray', trxAntenna);
            obj.txPower = txPower;
            % initialize precoder
            if isa(trxAntenna, "phased.ULA")
                obj.precoding = ones(1, trxAntenna.NumElements);
            elseif isa(trxAntenna, "phased.URA")
                obj.precoding = ones(trxAntenna.Size);
            else
                error('Unsupported antenna type');
            end
        end

        function setPrecoding(obj, fc, angle)
            obj.precoding = obj.steeringObject(fc, angle);
        end
    end
end