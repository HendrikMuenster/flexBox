classdef HuberProxDual < basicProx
    properties
    end
    
    methods
        function obj = HuberProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                %factor = 1 ./ (1+main.params.sigma{dualNumbers(i)}*tmpFactor);
                
                %main.yTilde{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ (1+main.params.sigma{dualNumbers(i)}*obj.epsi / obj.factor);
                main.yTilde{dualNumbers(i)} = main.yTilde{dualNumbers(i)} .* obj.factor ./ (obj.factor + obj.epsi.*main.params.sigma{dualNumbers(i)});
				norm = norm + (main.yTilde{dualNumbers(i)}).^2;
            end
            norm = max(1,sqrt(norm) / obj.factor);
            
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ norm;
            end
        end
        
    end
end