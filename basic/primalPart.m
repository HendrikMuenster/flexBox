classdef primalPart < functionalPart
    
    properties(SetAccess = protected, GetAccess = public)
    end
    
    methods
        function obj = primalPart(alpha)
            obj = obj@functionalPart(alpha);
        end
        
        function init(obj)
            
        end
    end
    
    methods (Abstract)
        applyProx(obj)
    end
end