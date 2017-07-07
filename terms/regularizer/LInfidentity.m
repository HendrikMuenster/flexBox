classdef LInfidentity < basicIdentity & LInfProxDual
    properties
    end

    methods
        function obj = LInfidentity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
        end

    end
end
