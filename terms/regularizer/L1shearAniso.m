%prox for F = alpha / 2 |\nabla u|^2
classdef L1shearAniso < basicShear & L1AnisoProxDual
    properties
    end
    
    methods
        function obj = L1shearAniso(alpha,dims,varargin)
            obj = obj@basicShear(alpha,dims,varargin);
        end
    end
end