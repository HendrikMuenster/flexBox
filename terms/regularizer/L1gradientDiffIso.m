%represents the term
%\alpha |\nabla(u-v)|_1,2
%in the L_1,2 norm
%corresponds to two primal variables (u,v)
classdef L1gradientDiffIso < basicGradientDifference & L1IsoProxDual
    properties
    end

    methods
        function obj = L1gradientDiffIso(alpha,dims,varargin)
            obj = obj@basicGradientDifference(alpha,dims,varargin);
        end
    end
end
