% 
classdef KLDataProxDual < handle
    properties
    end
    
    methods
        function obj = KLDataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = min(obj.factor,0.5*(obj.factor + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}+obj.factor).^2 + 4*(main.params.sigma{dualNumber}*obj.factor*obj.f{1} - obj.factor*main.yTilde{dualNumber}) )));
        end
        
    end
end