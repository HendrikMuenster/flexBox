% prox for point-wise projection onto L2 balls
% min_y 1/2\|y-x\|_2^2 + \delta_{ \|y/\alpha + b\|_2 \leq 1 }
classdef L1IsoProxDualShift < basicProx
    properties
    end
    
    methods
        function obj = L1IsoProxDualShift()
        end

        function applyProx(obj,main,dualNumbers,~)
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                norm = norm + (main.yTilde{dualNumbers(i)}/obj.factor - obj.b{i}).^2;
            end
            norm = max(1,sqrt(norm));
            
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = obj.factor *(main.yTilde{dualNumbers(i)}/obj.factor - obj.b{i}) ./ norm + obj.factor*obj.r{i};
            end
        end
        
    end
end