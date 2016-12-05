%represents the term
%\alpha |K(u,v)|_1,2 where K(u,v) = [u_x + v_y;u_x - v_y]
%in the L_1,2 norm defined as \sum abs( sqrt( (u_x + v_y)^2 + (u_x - v_y)^2 ) )
%corresponds to two primal variables (u,v)
classdef L1divAnticurlIso < basicDivAnticurl & L1IsoProxDual
    properties
    end

    methods
        function obj = L1divAnticurlIso(alpha,dims,varargin)
            obj = obj@basicDivAnticurl(alpha,dims,varargin);
        end
    end
end
