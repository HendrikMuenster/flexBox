%
classdef L2gradient < basicGradient & L2proxDual
    properties
    end
    
    methods
        function obj = L2gradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
			
			obj.CPPsupport = 1;
        end
    end
end