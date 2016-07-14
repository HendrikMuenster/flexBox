%
classdef huberOperator < basicDualizedOperator & L1HuberProxDual
    properties
        epsi
    end
    
    methods
        function obj = huberOperator(alpha,numPrimals,A,epsi,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
            obj.epsi = epsi;
        end
        
    end
end