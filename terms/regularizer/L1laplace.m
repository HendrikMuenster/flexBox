classdef L1laplace < basicLaplace & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1laplace(alpha,dims,varargin)
            obj = obj@basicLaplace(alpha,dims,varargin);
        end
    end
end