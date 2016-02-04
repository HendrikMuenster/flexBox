%prox for F = alpha / 2 |\nabla u|^2
classdef L1hessian < basicHessian & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1hessian(alpha,dims,varargin)
            obj = obj@basicHessian(alpha,dims,varargin);
        end
    end
end