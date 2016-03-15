%
classdef L1divergence < basicDivergence & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1divergence(alpha,dims,varargin)
            obj = obj@basicDivergence(alpha,dims,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end