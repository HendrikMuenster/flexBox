%prox for F = alpha / 2 |\nabla u|^2
classdef L1shearIso < basicShear & L1IsoProxDual
    properties
    end
    
    methods
        function obj = L1shearIso(alpha,dims,varargin)
            obj = obj@basicShear(alpha,dims,varargin);
        end
    end
end