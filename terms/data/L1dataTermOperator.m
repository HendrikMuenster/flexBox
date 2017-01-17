%term for \alpha |Au-f|_1
classdef L1dataTermOperator < basicDualizedDataterm & L1DataProxDual
    methods
        function obj = L1dataTermOperator(alpha,A,f,varargin)
            if (iscell(A))
                numPrimals = numel(A);
            else
                numPrimals = 1;
            end
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,f,varargin);
        end
    end
end