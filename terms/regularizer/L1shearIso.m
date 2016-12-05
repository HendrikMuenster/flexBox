%represents the term
%\alpha |K(u,v)|_1,2 where K(u,v) =  is the shear operator
%in the L_1,2 norm defined as \sum abs( abs(  ) )
%corresponds to two primal variables (u,v)
classdef L1shearIso < basicShear & L1IsoProxDual
    properties
    end

    methods
        function obj = L1shearIso(alpha,dims,varargin)
            obj = obj@basicShear(alpha,dims,varargin);
        end
    end
end
