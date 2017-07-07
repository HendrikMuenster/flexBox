classdef LInfgradient < basicGradient & LInfProxDual
    properties
    end

    methods
        function obj = LInfgradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end
