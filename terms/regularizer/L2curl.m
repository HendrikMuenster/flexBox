%represents the term
%\alpha/2 |K(u,v)|_2^2 where K(u,v) = [u_y - v_x] is the curl operator
%in the L_2^2 norm defined as \sum ( u_y - v_x ).^2
%corresponds to two primal variables (u,v)
classdef L2curl < basicCurl & L2proxDual
    properties
    end

    methods
        function obj = L2curl(alpha,dims,varargin)
            obj = obj@basicCurl(alpha,dims,varargin);
        end

    end
end
