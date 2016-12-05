%represents the term
%\alpha |K(u,v)|_2^2 where K(u,v) =  is the shear operator
%in the L_2^2 norm defined as \sum ().^2 
%corresponds to two primal variables (u,v)
classdef L2shear < basicShear & L2proxDual
    properties
    end

    methods
        function obj = L2shear(alpha,dims,varargin)
            obj = obj@basicShear(alpha,dims,varargin);
        end
    end
end
