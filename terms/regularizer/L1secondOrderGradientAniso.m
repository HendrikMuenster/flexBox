%represents the term
%\alpha |\nabla u - v|_1,1
%in the L_1,1 norm
%corresponds to two primal variables (u,v)
classdef L1secondOrderGradientAniso < basicSecondOrderGradient & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1secondOrderGradientAniso(alpha,dims,varargin)
            obj = obj@basicSecondOrderGradient(alpha,dims,varargin);
        end
    end
end
