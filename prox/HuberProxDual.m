classdef HuberProxDual < handle
    properties
    end
    
    methods
        function obj = HuberProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            
            tmpFactor = obj.epsi / obj.factor;
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                %factor = 1 ./ (1+main.params.sigma{dualNumbers(i)}*tmpFactor);
                
                main.yTilde{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ (1+main.params.sigma{dualNumbers(i)}*tmpFactor);
                norm = norm + (main.yTilde{dualNumbers(i)}).^2;
            end
            norm = max(1,sqrt(norm) / obj.factor);
            
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ norm;
            end
        end
        
    end
end