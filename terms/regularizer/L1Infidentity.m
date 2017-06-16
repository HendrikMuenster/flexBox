classdef L1Infidentity < basicIdentity & L1InfProxDual
    properties
    end

    methods
        function obj = L1Infidentity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
        end

    end
end
