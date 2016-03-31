%prox for G = alpha  |u-f|
classdef L1dataTermOld < primalPart
    
    properties
        f
    end

    methods
        function obj = L1dataTermOld(alpha,f)
            obj = obj@primalPart(alpha);%only one primal variable
            obj.f = f(:);
			
			obj.CPPsupport = 1;
        end
           
        function applyProx(obj,main,primalNumber)
            lambda = main.params.tau{primalNumber} * obj.factor;
            
            c1 = (main.xTilde{primalNumber} - obj.f) < - lambda;
            c2 = (main.xTilde{primalNumber} - obj.f) > lambda;
            
            main.x{primalNumber} = main.xTilde{primalNumber} + lambda*(c1 - c2) - (main.xTilde{primalNumber} - obj.f).*(1-c1-c2);
        end
    end
end