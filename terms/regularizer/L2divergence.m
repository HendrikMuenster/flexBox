%represents the term
%\alpha |K(u,v)|_2^2 where K(u,v) = [u_x + v_y]
%in the L_2^2 norm defined as \sum  ( u_x + v_y ).^2
%corresponds to two primal variables (u,v)
classdef L2divergence < basicDivergence & L2proxDual
    properties
    end

    methods
        function obj = L2divergence(alpha,dims,varargin)
            obj = obj@basicDivergence(alpha,dims,varargin);
        end
    end
end
