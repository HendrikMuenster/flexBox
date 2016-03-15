%
classdef L1operatorAniso < basicDualizedOperator & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1operatorAniso(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end