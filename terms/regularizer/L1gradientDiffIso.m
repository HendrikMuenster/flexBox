%
classdef L1gradientDiffIso < basicGradientDifference & L1IsoProxDual
    properties
    end
    
    methods
        function obj = L1gradientDiffIso(alpha,dims,varargin)
            obj = obj@basicGradientDifference(alpha,dims,varargin);
			
			obj.CPPsupport = 1;
        end
    end
end