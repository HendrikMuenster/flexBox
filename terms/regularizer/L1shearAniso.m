%represents the term
%\alpha |K(u,v)|_1,1 where K(u,v) =  is the shear operator
%in the L_1,1 norm defined as \sum abs( abs(  ) )
%corresponds to two primal variables (u,v)
classdef L1shearAniso < basicShear & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1shearAniso(alpha,dims,varargin)
            obj = obj@basicShear(alpha,dims,varargin);
        end
    end
end
