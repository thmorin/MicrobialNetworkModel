classdef parameters
    properties 
        Ks
        Kd
        mu_max
        Y
        kN
        kP
    end
    methods
        function params = parameters(K,death,max,yield,kN,kP)
            if nargin == 6
                params.Ks     = K;
                params.Kd     = death;
                params.mu_max = max;
                params.Y      = yield;
                params.kN     = kN;
                params.kP     = kP;
            end
        end
    end
end

                