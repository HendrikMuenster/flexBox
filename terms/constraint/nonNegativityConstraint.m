%
classdef nonNegativityConstraint < basicDualizedDataterm & constraintBoxDualized
    methods
        function obj = nonNegativityConstraint(sizeVariable)
            obj = obj@basicDualizedDataterm(1,1,identityOperator(prod(sizeVariable)),cell(1,1));
            
            obj.minVal = 0;
            obj.maxVal = 1e10;
        end
    end
end