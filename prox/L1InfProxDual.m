classdef L1InfProxDual < basicProx
    properties
    end

    methods
        function obj = L1InfProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            for i=1:obj.numVars
                currentY = main.yTilde{dualNumbers(i)};

                %project component-wise onto l1 ball [Author: John Duchi (jduchi@cs.berkeley.edu)]
                if (norm(currentY, 1) < obj.factor)
                  projectedY = currentY;
                else
                  u = sort(abs(currentY), 'descend');
                  sv = cumsum(u);
                  rho = find(u > (sv - obj.factor) ./ (1:length(u)), 1, 'last');
                  theta = max(0, (sv(rho) - obj.factor) / rho);
                  w = sign(currentY) .* max(abs(currentY) - theta, 0);
                end

                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} - obj.factor * projectedY;
            end
        end
    end
end
