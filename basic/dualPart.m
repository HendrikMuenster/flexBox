classdef dualPart < functionalPart
    properties
        operator  %operator
        operatorT %transposed operator
        numVars   %num of dual variables
        length    %vector of length numVars containing the length of corresponding dual variable
        myTau     %primal dual normalization factor for primal part
        mySigma   %primal dual normalization factor for dual part
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
