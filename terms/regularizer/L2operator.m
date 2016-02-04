%prox for F = alpha / 2 |\nabla u|^2
classdef L2operator < basicDualizedOperator & L2proxDual
    properties
    end
    
    methods
        function obj = L2operator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end