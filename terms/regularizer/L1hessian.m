%represents the term
%\alpha |Hu|_1,1 where H is the hessian operator
%in the L_1,1 norm
%corresponds to one primal variable u
classdef L1hessian < basicHessian & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1hessian(alpha,dims,varargin)
            obj = obj@basicHessian(alpha,dims,varargin);
        end
    end
end
