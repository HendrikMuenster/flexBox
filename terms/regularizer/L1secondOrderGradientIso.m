%represents the term
%\alpha |\nabla u - v|_1,2
%in the L_1,2 norm
%corresponds to two primal variables (u,v)
classdef L1secondOrderGradientIso < basicSecondOrderGradient & L1IsoProxDual
    properties
    end

    methods
        function obj = L1secondOrderGradientIso(alpha,dims,varargin)
            obj = obj@basicSecondOrderGradient(alpha,dims,varargin);
        end
    end
end
