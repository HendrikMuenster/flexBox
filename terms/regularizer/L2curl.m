%corresponds to two primal variables (u,v) and represents the curl operator 
%K(u,v) = [u_y - v_x] in the L_2^2 norm defined as
% \sum ( u_y - v_x ).^2
classdef L2curl < basicCurl & L2proxDual
    properties
    end
    
    methods
        function obj = L2curl(alpha,dims,varargin)
            obj = obj@basicCurl(alpha,dims,varargin);
            
            obj.CPPsupport = 1;
        end

    end
end