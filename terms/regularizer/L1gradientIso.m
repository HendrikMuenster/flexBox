%
classdef L1gradientIso < basicGradient & L1IsoProxDual
    properties
    end
    
    methods
        function obj = L1gradientIso(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
			
			obj.CPPsupport = 1;
        end
    end
end