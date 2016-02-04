%corresponds to two primal variables (u,v) and prepresents the operator 
%K(u,v) = [u_x + v_y;u_y - v_x] in the L_1,2 norm defined as
% \sum abs( sqrt( (u_x + v_y)^2 + (u_y - v_x)^2 ) )
classdef L1divCurlIso < basicDivCurl & L1IsoProxDual
    properties
    end
    
    methods
        function obj = L1divCurlIso(alpha,dims,varargin)
            obj = obj@basicDivCurl(alpha,dims,varargin);
        end
    end
end