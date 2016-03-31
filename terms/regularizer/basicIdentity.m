classdef basicIdentity < basicDualizedOperator
    properties
    end
    
    methods
        function obj = basicIdentity(alpha,dims,varargin)
            obj = obj@basicDualizedOperator(alpha,1,identityOperator(prod(dims)),varargin);
        end
    end
end