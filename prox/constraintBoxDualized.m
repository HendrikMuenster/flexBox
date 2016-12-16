classdef constraintBoxDualized < handle
    properties
        minVal
        maxVal
    end
    
    methods
        function obj = constraintBoxDualized()
            %default values, change in class if necessary
            obj.minVal = 0;
            obj.maxVal = 1;
        end

        function applyProx(obj,main,dualNumbers,~)
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = max(0,main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)} * obj.maxVal) + min(0,main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)} * obj.minVal);
                %main.yTilde{dualNumbers(i)} .* obj.maxVal .* (main.yTilde{dualNumbers(i)} > 0) + main.yTilde{dualNumbers(i)} .* obj.minVal .* (main.yTilde{dualNumbers(i)} < 0);
            end
        end
        
    end
end