% 
classdef L1HuberDataProxDual < handle
    properties
    end
    
    methods
        function obj = L1HuberDataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            factor = 1 ./ (1+main.params.sigma{dualNumber}*obj.epsi ./ obj.factor);

            main.y{dualNumber} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber}.* factor - factor*main.params.sigma{dualNumber}*obj.f(:)));
        end
        
    end
end