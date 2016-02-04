classdef MinMaxProxPrimal < handle
    properties
        minVal
        maxVal
    end
    
    methods
        function obj = MinMaxProxPrimal()
            %default values, change in class if necessary
            obj.minVal = 0;
            obj.maxVal = 1;
        end

        function applyProx(obj,main,primalNumbers)
            for i=1:obj.numVars
                main.x{primalNumbers(i)} = min(obj.maxVal,max(obj.minVal,main.xTilde{primalNumbers(i)}));
            end
        end
        
    end
end