%prox for F = alpha / 2 |\nabla u|^2
classdef L2curl < basicCurl & L2proxDual
    properties
    end
    
    methods
        function obj = L2curl(alpha,dims,varargin)
            obj = obj@basicCurl(alpha,dims,varargin);
        end
        
        %function call to preallocate usefull things
        function init(varargin)
        end

    end
end