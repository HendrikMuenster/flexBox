%represents the term
%\alpha |K(u,v)|_1,1 where K(u,v) = [u_x^T - v_y^T] is the switched curl operator
%in the L_1,1 norm defined as \sum abs( abs( u_x^T - v_y^T ) )
%corresponds to two primal variables (u,v)
classdef L1curlSwitched < basicCurl & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1curlSwitched(alpha,dims,varargin)
            obj = obj@basicCurl(alpha,dims,'switchedCurl',1);
        end
    end
end
