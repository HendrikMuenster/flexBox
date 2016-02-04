classdef functionalPart < handle
    
    properties
        factor
		CPPsupport;
        numPrimals %corresponding number of primals
    end
    
    methods
        function obj = functionalPart(alpha)
            obj.factor = alpha;
			obj.CPPsupport = 0;
            obj.numPrimals = 1;
        end
        
        %function call to preallocate usefull things
        function init(varargin)
        end
    end
    

    
    methods (Abstract)
        applyProx(obj)
    end
end