%
classdef boxConstraint < basicDualizedDataterm & constraintBoxDualized
    methods
        function obj = boxConstraint(u1,u2,dims)
            obj = obj@basicDualizedDataterm(1,1,identityOperator(prod(dims)),cell(1,1));
            
            obj.minVal = u1;
            obj.maxVal = u2;
        end
    end
end