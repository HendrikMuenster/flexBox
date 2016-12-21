%represents the term
%\alpha |\nabla u|_1,1
%in the L_1,1 norm
%corresponds to one primal variable u
classdef L1gradientAniso < basicGradient & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1gradientAniso(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end
