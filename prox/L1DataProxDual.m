% prox for point-wise projection onto L1 balls
% min_y 1/2\|y-x\|_2^2 + <y,f> + \delta_{ \|y\|_2 \leq \alpha }
classdef L1DataProxDual < basicProx
    properties
    end
    
    methods
        function obj = L1DataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            for i=1:obj.numVars
                main.y{dualNumber(i)} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber(i)} - main.params.sigma{dualNumber(i)}.*obj.f{i}(:)));
            end
        end
    end
end