% prox for Frobenius norm
% min_y 1/2\|y-x\|_2^2 + \delta_{ \|y\|_F \leq \alpha }
classdef FrobeniusProxDual < handle
    properties
    end
    
    methods
        function obj = FrobeniusProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                norm = norm + sum(main.yTilde{dualNumbers(i)}.^2);
            end
            norm = max(1,sqrt(norm)./obj.factor);
            
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ norm;
            end
        end
        
    end
end