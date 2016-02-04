%
classdef L1gradientAniso < basicGradient & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1gradientAniso(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
            
			obj.CPPsupport = 1;
        end
    end
end