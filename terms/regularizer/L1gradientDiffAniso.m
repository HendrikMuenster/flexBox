%
classdef L1gradientDiffAniso < basicGradientDifference & L1anisoProxDual
    properties
    end
    
    methods
        function obj = L1gradientDiffAniso(alpha,dims,varargin)
            obj = obj@basicGradientDifference(alpha,dims,varargin);
			
			obj.CPPsupport = 1;
        end
    end
end