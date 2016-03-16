%corresponds to two primal variables (u,v) and represents the divergence operator 
%K(u,v) = [u_x + v_y] in the L_1 norm defined as
% \sum abs( u_x + v_y )
classdef L1divergence < basicDivergence & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1divergence(alpha,dims,varargin)
            obj = obj@basicDivergence(alpha,dims,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end