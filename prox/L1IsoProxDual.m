% prox for point-wise projection onto L2 balls
% min_y 1/2\|y-x\|_2^2 + \delta_{ \|y\|_2 \leq \alpha }
classdef L1IsoProxDual < handle
    properties
    end
    
    methods
        function obj = L1IsoProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                norm = norm + main.yTilde{dualNumbers(i)}.^2;
            end
            norm = max(1,sqrt(norm)./obj.factor);
            
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ norm;
            end
        end
        
    end
end