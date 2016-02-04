%
classdef L2identity < basicIdentity & L2proxDual
    properties
    end
    
    methods
        function obj = L2identity(alpha,dims)
            obj = obj@ basicIdentity(alpha,dims);
        end

        %function call to preallocate usefull things
        function init(obj,myNumber,main)
        end
        
    end
end