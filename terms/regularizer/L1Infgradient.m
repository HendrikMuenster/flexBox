classdef L1Infgradient < basicGradient & L1InfProxDual
    properties
    end

    methods
        function obj = L1Infgradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end
