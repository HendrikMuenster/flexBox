% prox for L2 with additional f
% min_y 1/2\|y-x\|_2^2 + <y,f> + \alpha\|y\|_2^2 
classdef L2DataProxDual < handle
    properties
    end
    
    methods
        function obj = L2DataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = (obj.factor/(main.params.sigma{dualNumber}+obj.factor)) * (main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:));
        end
    end
end