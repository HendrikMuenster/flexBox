%
classdef L1secondOrderGradientAniso < basicSecondOrderGradient & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1secondOrderGradientAniso(alpha,dims,varargin)
            obj = obj@basicSecondOrderGradient(alpha,dims,varargin);
            obj.CPPsupport = 1;
        end
    end
end