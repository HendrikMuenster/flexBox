%basic abstract operator class
classdef basicOperator

    properties
    end
    
    methods (Abstract)
        mtimes(obj,vector)
        returnMatrix(obj)
    end
    
end

