classdef dualPart < functionalPart
    
    properties
        operator
        numVars
        length
        myTau
        mySigma
    end
    
    methods
        function obj = dualPart(alpha)
            obj = obj@functionalPart(alpha);
        end
    end
    
    methods (Abstract)
        applyProx(obj)
    end
end