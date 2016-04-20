%
classdef huberGradient < basicGradient & L1HuberProxDual
    properties
        epsi
    end
    
    methods
        function obj = huberGradient(alpha,dims,epsi,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
            obj.epsi = epsi;
        end
        
    end
end