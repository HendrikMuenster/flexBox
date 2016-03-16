%corresponds to two primal variables (u,v) and represents the divergence operator 
%K(u,v) = [u_x + v_y] in the L_2^2 norm defined as
% \sum  ( u_x + v_y ).^2 
classdef L2divergence < basicDivergence & L2proxDual
    properties
    end
    
    methods
        function obj = L2divergence(alpha,dims,varargin)
            obj = obj@basicDivergence(alpha,dims,varargin);
            
            obj.CPPsupport = 1;
        end
        
        %function call to preallocate usefull things
        function init(varargin)
        end
    end
end