classdef L2Infoperator < basicDualizedOperator & L2InfProxDual
    properties
    end

    methods
        function obj = L2Infoperator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
