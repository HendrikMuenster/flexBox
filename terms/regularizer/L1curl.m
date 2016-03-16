%corresponds to two primal variables (u,v) and represents the curl operator 
%K(u,v) = [u_y^T - v_x^T] in the L_1 norm defined as
% \sum abs( abs( u_y^T - v_x^T ) )
classdef L1curl < basicCurl & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1curl(alpha,dims,varargin)
            obj = obj@basicCurl(alpha,dims,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end