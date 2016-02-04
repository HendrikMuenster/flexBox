classdef primalVar < handle
    properties
        dims
        var
        old
    end
    
    methods
        function obj = primalVar(dims)
            obj.var = zeros(prod(dims),1);
            obj.old = obj.var;
            obj.dims = dims;
        end
        
        function result = get(obj)
            result = obj.var;
        end
        
        function result = getBar(obj)
            result = 2*obj.var - 2*obj.old;
        end
        
        function set(obj,value)
            obj.old = obj.var;
            obj.var = value;
        end
    end
    
end

