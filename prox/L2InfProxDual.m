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
            sortyTildeNorm = unique(sort(yTildeNorm));
            
            h = @(lambda) (yTildeNorm > lambda) .* (yTildeNorm - lambda);
            g = @(lambda) sum(h(lambda)) - obj.factor;
            lambda = max(0,fzero(g, 0));
            
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = (yTildeNorm > lambda) .* main.yTilde{dualNumbers(i)} .* (1 - lambda ./ yTildeNorm);
                main.y{dualNumbers(i)}(yTildeNorm <= lambda) = 0;
            end
            
            

%             %first guess for main.y
%             for i=1:obj.numVars
%                 main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)};
%             end
%             yNorm = yTildeNorm;
%             
%             counter = 1;
%             while (sum(yNorm) >= obj.factor) && (counter <= length(sortyTildeNorm))
%                 lambda = sortyTildeNorm(counter);
%                 
%                 %update guess for main.y
%                 for i=1:obj.numVars
%                     main.y{dualNumbers(i)} = (yTildeNorm > lambda) .* main.yTilde{dualNumbers(i)} .* (1 - lambda ./ yTildeNorm);
%                     main.y{dualNumbers(i)}(yTildeNorm <= lambda) = 0;
%                 end
%                 
%                 yNorm = 0;
%                 for i=1:obj.numVars
%                     yNorm = yNorm + main.y{dualNumbers(i)}.^2;
%                 end
%                 yNorm = sqrt(yNorm);
%                 %yNorm = yTildeNorm - lambda;
%                 
%                 counter = counter + 1;
%             end
        end       
    end
end
