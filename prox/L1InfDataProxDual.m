classdef L1InfDataProxDual < basicProx
    properties
    end

    methods
        function obj = L1InfDataProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                norm = norm + abs(main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)}.*obj.f{i}(:));
            end
            norm = max(obj.factor,norm);

            for i=1:obj.numVars
                main.y{dualNumbers(i)} = obj.factor * (main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)}.*obj.f{i}(:)) ./ norm;
            end
        end
    end
end
