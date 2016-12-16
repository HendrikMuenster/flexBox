%
classdef nonPositivityConstraint < basicDualizedDataterm & constraintBoxDualized
    methods
        function obj = nonPositivityConstraint(sizeVariable)
            obj = obj@basicDualizedDataterm(1,1,identityOperator(prod(sizeVariable)),cell(1,1));
            
            obj.minVal = -1e10;
            obj.maxVal = 0;
        end
    end
end