classdef LInfdataTermOperator < basicDualizedDataterm & LInfDataProxDual
    methods
        function obj = LInfdataTermOperator(alpha,A,f,varargin)
            if (iscell(A))
                numPrimals = numel(A);
            else
                numPrimals = 1;
            end
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,f,varargin);
        end
    end
end
