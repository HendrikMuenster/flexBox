%
classdef L2divergence < basicDivergence & L2proxDual
    properties
    end
    
    methods
        function obj = L2divergence(alpha,dims,varargin)
            obj = obj@basicDivergence(alpha,dims,varargin);
        end
        
        %function call to preallocate usefull things
        function init(varargin)
        end
    end
end