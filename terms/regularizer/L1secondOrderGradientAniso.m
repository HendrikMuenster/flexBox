%
classdef L1secondOrderGradientAniso < basicSecondOrderGradient & L1anisoProxDual
    properties
    end
    
    methods
        function obj = L1secondOrderGradientAniso(alpha,dims,varargin)
            obj = obj@basicSecondOrderGradient(alpha,dims,varargin);
        end
    end
end