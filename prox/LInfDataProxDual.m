classdef LInfDataProxDual < basicProx
    properties
    end

    methods
        function obj = LInfDataProxDual()
    end

    function applyProx(obj,main,dualNumbers,~)
        for i=1:obj.numVars
            %project main.yTilde{dualNumbers(i)} onto l1 balls of radius alpha
            v = main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)} .* obj.f{i};
            b = obj.factor;
            if (norm(v, 1) < b)
                main.y{dualNumbers(i)} = v;
            else
                u = sort(abs(v),'descend');
                sv = cumsum(u);
                rho = find(u > (sv - b) ./ ((1:length(u))'), 1, 'last');
                theta = max(0, (sv(rho) - b) / rho);
                main.y{dualNumbers(i)} = sign(v) .* max(abs(v) - theta, 0);
            end
        end
    end
    end
end
