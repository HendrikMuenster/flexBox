%represents the term
%\alpha |Hu|_2^2 where H is the hessian operator
%in the L_2^2 norm
%corresponds to one primal variable u
classdef L2hessian < basicHessian & L2proxDual
    properties
    end

    methods
        function obj = L2hessian(alpha,dims,varargin)
            obj = obj@basicHessian(alpha,dims,varargin);
        end
    end
end
