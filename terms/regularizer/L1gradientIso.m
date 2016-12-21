%represents the term
%\alpha |\nabla u|_1,2
%in the L_1,2 norm
%corresponds to one primal variable u
classdef L1gradientIso < basicGradient & L1IsoProxDual
    properties
    end

    methods
        function obj = L1gradientIso(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end
