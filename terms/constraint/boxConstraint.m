%
classdef boxConstraint < basicDualizedDataterm & constraintBoxDualized
    methods
        function obj = boxConstraint(minVal,maxVal,sizeVariable)
            obj = obj@basicDualizedDataterm(1,1,identityOperator(prod(sizeVariable)),cell(1,1));
            
            obj.minVal = minVal;
            obj.maxVal = maxVal;
        end
    end
end