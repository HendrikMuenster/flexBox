%
classdef L2laplace < basicLaplace & L2proxDual
    properties
    end
    
    methods
        function obj = L2laplace(alpha,dims,varargin)
            obj = obj@basicLaplace(alpha,dims,varargin);
        end
    end
end