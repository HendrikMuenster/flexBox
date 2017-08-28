classdef L2Infgradient < basicGradient & L2InfProxDual
    properties
    end

    methods
        function obj = L2Infgradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end
