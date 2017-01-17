%term for \alpha/2 |Au-f|_2^2
classdef L2dataTermOperator < basicDualizedDataterm & L2DataProxDual
    methods
        function obj = L2dataTermOperator(alpha,A,f,varargin)
            if (iscell(A))
                numPrimals = numel(A);
            else
                numPrimals = 1;
            end
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,f,varargin);
        end
    end
end