%corresponds to two primal variables (u,v) and represents the switched curl
%operator operator K(u,v) = [u_x^T - v_y^T] in the L_1,1 norm defined as
% \sum abs( abs( u_x^T - v_y^T ) )
classdef L1curlSwitched < basicCurl & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1curlSwitched(alpha,dims,varargin)
            obj = obj@basicCurl(alpha,dims,'switchedCurl',1);
        end
    end
end