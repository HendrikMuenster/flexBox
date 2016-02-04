%
classdef L1gradientDifferenceIso < basicGradientDifference & L1IsoProxDual
    properties
    end
    
    methods
        function obj = L1gradientDifferenceIso(alpha,dims,varargin)
            obj = obj@basicGradientDifference(alpha,dims,varargin);
			
			obj.CPPsupport = 1;
        end
    end
end