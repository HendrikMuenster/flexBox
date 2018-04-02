classdef nuclearIdentity < basicIdentity & nuclearProxDual
    properties
    end

    methods
        function obj = nuclearIdentity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
        end

    end
end
