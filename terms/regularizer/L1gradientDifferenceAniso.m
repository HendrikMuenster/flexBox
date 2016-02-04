%
classdef L1gradientDifferenceAniso < basicGradientDifference & L1anisoProxDual
    properties
    end
    
    methods
        function obj = L1gradientDifferenceAniso(alpha,dims,varargin)
            obj = obj@basicGradientDifference(alpha,dims,varargin);
			
			obj.CPPsupport = 1;
        end
    end
end