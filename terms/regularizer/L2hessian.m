%prox for F = alpha / 2 |\nabla u|^2
classdef L2hessian < basicHessian & L2proxDual
    properties
    end
    
    methods
        function obj = L2hessian(alpha,dims,varargin)
            obj = obj@basicHessian(alpha,dims,varargin);
        end
    end
end