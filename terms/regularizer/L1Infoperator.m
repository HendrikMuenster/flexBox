classdef L1Infoperator < basicDualizedOperator & L1InfProxDual
    properties
    end

    methods
        function obj = L1Infoperator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
