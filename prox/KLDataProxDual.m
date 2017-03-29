% 
classdef KLDataProxDual < basicProx
    properties
    end
    
    methods
        function obj = KLDataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            
            
            %solution1 = 0.5*(obj.factor + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}+obj.factor).^2 + 4*(main.params.sigma{dualNumber}*obj.factor*obj.f{1} - obj.factor*main.yTilde{dualNumber}) ));
            %solution2 = 0.5*(obj.factor + main.yTilde{dualNumber} + sqrt( (main.yTilde{dualNumber}+obj.factor).^2 + 4*(main.params.sigma{dualNumber}*obj.factor*obj.f{1} - obj.factor*main.yTilde{dualNumber}) ));
            

            main.y{dualNumber} = 0.5*(obj.factor + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}+obj.factor).^2 + 4*(main.params.sigma{dualNumber}.*obj.factor.*obj.f{1} - obj.factor*main.yTilde{dualNumber}) ));
        end
        
    end
end