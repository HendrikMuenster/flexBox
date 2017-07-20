% 
classdef HuberDataProxDual < basicProx
    properties
    end
    
    methods
        function obj = HuberDataProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            %calc norm
            norm = 0;
            for i=1:obj.numVars
                
				main.yTilde{dualNumbers(i)} = (main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)}.*obj.f{i}(:)) .* obj.factor ./ (obj.factor + obj.epsi.*main.params.sigma{dualNumbers(i)});
				
				norm = max(1,abs(main.yTilde{dualNumbers(i)}) / obj.factor);
				
				main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} ./ norm;
            end

        end
		
		
		
        
    end
end