%represents the term
%\alpha |\nabla(u-v)|_1,1
%in the L_1,1 norm
%corresponds to two primal variables (u,v)
classdef L1gradientDiffAniso < basicGradientDifference & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1gradientDiffAniso(alpha,dims,varargin)
            obj = obj@basicGradientDifference(alpha,dims,varargin);
        end
    end
end
