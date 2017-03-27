classdef IdentityProxPrimal < basicProx
    properties
    end
    
    methods
        function obj = IdentityProxPrimal()
        end

        function applyProx(obj,main,primalNumbers)
            for i=1:numel(primalNumbers)
                main.x{primalNumbers(i)} = main.xTilde{primalNumbers(i)};
            end
        end
        
    end
end