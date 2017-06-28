classdef L1InfDataProxDual < basicProx
    properties
    end

    methods
        function obj = L1InfDataProxDual()
    end

    function applyProx(obj,main,dualNumbers,~)
        for i=1:obj.numVars
            currentY = main.yTilde{dualNumbers(i)} ./ main.params.sigma{dualNumbers(i)} ...
                - main.params.sigma{dualNumbers(i)} .* obj.f{i};


            %project component-wise onto l1 ball [Author: John Duchi (jduchi@cs.berkeley.edu)]
            if (norm(currentY, 1) < obj.factor)
              projectedY = currentY;
            else
              u = sort(abs(currentY), 'descend');
              sv = cumsum(u);
              rho = find(u > ((sv - obj.factor) ./ (1:length(u))'), 1, 'last');
              theta = max(0, (sv(rho) - obj.factor) / rho);
              projectedY = sign(currentY) .* max(abs(currentY) - theta, 0);
            end
            main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)} .* projectedY;
        end
    end
    end
end
