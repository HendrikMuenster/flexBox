% prox for point-wise projection onto L1 balls
% min_y 1/2\|y-x\|_2^2 + <y,f> + \delta_{ \|y\|_2 \leq \alpha }
classdef L1DataProxDual < handle
    properties
    end
    
    methods
        function obj = L1DataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:)));
        end
        
    end
end