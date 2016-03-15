classdef basicIdentity < basicDualizedOperator & tildeMultiOperatorMultiDual
    properties
    end
    
    methods
        function obj = basicIdentity(alpha,dims,varargin)
            obj = obj@basicDualizedOperator(alpha,1,speye(prod(dims)),varargin);
        end
    end
end