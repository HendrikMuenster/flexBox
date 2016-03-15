%
classdef L1identity < basicIdentity & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1identity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
            
            obj.CPPsupport = 1;
        end
        
    end
end