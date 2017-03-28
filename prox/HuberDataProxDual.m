% 
classdef HuberDataProxDual < basicProx
    properties
    end
    
    methods
        function obj = HuberDataProxDual()
        end

        function applyProx(obj,main,dualNumber,~)
            %factor = 1 ./ (1+main.params.sigma{dualNumber}*obj.epsi ./ obj.factor);
            
            for i=1:obj.numVars
                tmpVar = obj.factor / (obj.factor + main.params.sigma{dualNumber(i)}*obj.epsi) * (main.yTilde{dualNumber(i)} - main.params.sigma{dualNumber(i)}*obj.f{i});
                
                main.y{dualNumber(i)} = tmpVar ./ max(1,abs(tmpVar) / obj.factor);
                
                %main.y{dualNumber(i)} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber(i)}.* factor - factor*main.params.sigma{dualNumber(i)}*obj.f{i}));
            end
        end
        
    end
end