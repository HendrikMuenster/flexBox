classdef L2proxDual < basicProx
    properties
    end
    
    methods
        function obj = L2proxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = obj.factor./(main.params.sigma{dualNumbers(i)}+obj.factor) .* main.yTilde{dualNumbers(i)};
            end
        end
        
    end
end