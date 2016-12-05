%represents the term
%\alpha |\Delta u|_2^2
%in the L_2^2 norm
%corresponds to one primal variable u
classdef L2laplace < basicLaplace & L2proxDual
    properties
    end

    methods
        function obj = L2laplace(alpha,dims,varargin)
            obj = obj@basicLaplace(alpha,dims,varargin);
        end
    end
end
