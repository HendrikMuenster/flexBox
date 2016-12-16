% 
classdef KLDataProxDual < handle
    properties
    end
    
    methods
        function obj = KLDataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = min(1,0.5*(1 + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}-1).^2 + 4*main.params.sigma{dualNumber}*obj.f{1} )));
        end
        
    end
end