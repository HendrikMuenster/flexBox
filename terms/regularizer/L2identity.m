%
classdef L2identity < basicIdentity & L2proxDual
    properties
    end
    
    methods
        function obj = L2identity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
            
            obj.CPPsupport = 1;
        end
        
    end
end