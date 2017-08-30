classdef L2InfProxDual < basicProx
    properties
    end

    methods
        function obj = L2InfProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            yTildeNorm = 0;
            for i=1:obj.numVars
                yTildeNorm = yTildeNorm + main.yTilde{dualNumbers(i)}.^2;
            end
            yTildeNorm = sqrt(yTildeNorm);

            sortyTildeNorm = sort(yTildeNorm, 'descend');
        
            yTildeSum = 0;
            g2 = -obj.factor;
            lambda = 0;
            for index=2:length(sortyTildeNorm)
                lambda = sortyTildeNorm(index);
                yTildeSum = yTildeSum + sortyTildeNorm(index - 1);
                g2 = yTildeSum - (index - 1) * lambda - obj.factor;
                
                if g2 >= 0
                    break;
                end
            end
            
            if g2 < 0
                lambda = 0;
            else
                lambda = (yTildeSum - obj.factor) / (index - 1);
            end
  
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = (yTildeNorm > lambda) .* main.yTilde{dualNumbers(i)} .* (1 - lambda ./ yTildeNorm);
                main.y{dualNumbers(i)}(yTildeNorm <= lambda) = 0;
            end
        end
    end
end
