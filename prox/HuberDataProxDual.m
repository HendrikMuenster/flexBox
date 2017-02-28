% 
classdef HuberDataProxDual < handle
    properties
    end
    
    methods
        function obj = HuberDataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            factor = 1 ./ (1+main.params.sigma{dualNumber}*obj.epsi ./ obj.factor);
            
            for i=1:obj.numVars
                main.y{dualNumber(i)} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber(i)}.* factor - factor*main.params.sigma{dualNumber(i)}*obj.f{i}));
            end
        end
        
    end
end