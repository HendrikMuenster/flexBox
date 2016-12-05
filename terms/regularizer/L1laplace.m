%represents the term
%\alpha |\Delta u|_1,1
%in the L_1,1 norm
%corresponds to one primal variable u
classdef L1laplace < basicLaplace & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1laplace(alpha,dims,varargin)
            obj = obj@basicLaplace(alpha,dims,varargin);
        end
    end
end
