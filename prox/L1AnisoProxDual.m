classdef L1AnisoProxDual < basicProx
    properties
    end
    
    methods
        function obj = L1AnisoProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumbers(i)}));
            end
        end
        
    end
end