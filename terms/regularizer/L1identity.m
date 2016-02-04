%
classdef L1identity < basicIdentity & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1identity(alpha,dims)
            obj = obj@ basicIdentity(alpha,dims);
        end

        %function call to preallocate usefull things
        function init(obj,myNumber,main)
        end
        
    end
end