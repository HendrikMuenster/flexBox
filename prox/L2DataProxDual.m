% prox for L2 with additional f
% min_y 1/2\|y-x\|_2^2 + <y,f> + \alpha\|y\|_2^2 
classdef L2DataProxDual < basicProx
    properties
    end
    
    methods
        function obj = L2DataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            for i=1:obj.numVars
                main.y{dualNumber(i)} = (obj.factor./(main.params.sigma{dualNumber(i)}+obj.factor)) .* (main.yTilde{dualNumber(i)} - main.params.sigma{dualNumber(i)}.*obj.f{i}(:));
            end
        end
    end
end