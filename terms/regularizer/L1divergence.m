%represents the term
%\alpha |K(u,v)|_1,1 where K(u,v) = [u_x + v_y]
%in the L_1,1 norm defined as \sum abs( u_x + v_y )
%corresponds to two primal variables (u,v)
classdef L1divergence < basicDivergence & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1divergence(alpha,dims,varargin)
            obj = obj@basicDivergence(alpha,dims,varargin);
        end
    end
end
