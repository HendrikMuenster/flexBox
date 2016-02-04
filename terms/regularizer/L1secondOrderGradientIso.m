%
classdef L1secondOrderGradientIso < basicSecondOrderGradient & L1IsoProxDual
    properties
    end
    
    methods
        function obj = L1secondOrderGradientIso(alpha,dims,varargin)
            obj = obj@basicSecondOrderGradient(alpha,dims,varargin);
        end
    end
end