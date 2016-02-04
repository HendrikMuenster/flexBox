%prox for F = alpha / 2 |\nabla u|^2
classdef L2shear < basicShear & L2proxDual
    properties
    end
    
    methods
        function obj = L2shear(alpha,dims,varargin)
            obj = obj@basicShear(alpha,dims,varargin);
        end
    end
end