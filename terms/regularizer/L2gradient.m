%represents the term
%\alpha |\nabla u|_2^2
%in the L_2^2 norm
%corresponds to one primal variable u
classdef L2gradient < basicGradient & L2proxDual
    properties
    end

    methods
        function obj = L2gradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end
