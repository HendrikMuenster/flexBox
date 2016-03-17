%prox for F = alpha / 2 |\nabla u|^2
classdef frobeniusOperator < basicDualizedOperator & FrobeniusProxDual
    properties
    end
    
    methods
        function obj = frobeniusOperator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end