classdef L2InfProxDual < basicProx
    properties
    end

    methods
        function obj = L2InfProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            %moreau L1IsoProxDual
            norm = 0;
            for i=1:obj.numVars
                norm = norm + main.yTilde{dualNumbers(i)}.^2;
            end
            norm = max(obj.factor,sqrt(norm));

            for i=1:obj.numVars
                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)} .* (obj.factor*main.yTilde{dualNumbers(i)} ./ norm);
            end
        end
    end
end
