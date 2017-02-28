%
classdef nonNegativityConstraintOperator < basicDualizedDataterm & constraintBoxDualized
    methods
        function obj = nonNegativityConstraintOperator(A)
            if (iscell(A))
                numPrimals = numel(A);
            else
                numPrimals = 1;
            end
		
            obj = obj@basicDualizedDataterm(1,numPrimals,A,cell(1,1));
            
            obj.minVal = 0;
            obj.maxVal = 1e10;
        end
    end
end