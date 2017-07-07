classdef LInfoperator < basicDualizedOperator & LInfProxDual
    properties
    end

    methods
        function obj = LInfoperator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
