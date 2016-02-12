%
classdef frobeniusGradient < basicGradient & FrobeniusProxDual
    properties
    end
    
    methods
        function obj = frobeniusGradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
        end
    end
end