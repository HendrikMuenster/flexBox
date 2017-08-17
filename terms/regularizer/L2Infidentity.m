classdef L2Infidentity < basicIdentity & L2InfProxDual
    properties
    end

    methods
        function obj = L2Infidentity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
        end

    end
end
