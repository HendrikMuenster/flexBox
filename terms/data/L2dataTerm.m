% L2dataTerm(dims,alpha,f) represents an L2 data term 
% alpha/2 ||u-f||_2^2
classdef L2dataTerm < primalPart
    
    properties
        f
        xIndices
    end

    methods
        function obj = L2dataTerm(alpha,f)
            obj = obj@primalPart(alpha);
            obj.f = f(:);
			
			obj.CPPsupport = 1;
        end
        
        function init(obj,myNumber,main)

        end
           
        function applyProx(obj,main,primalNumber)
            main.x{primalNumber} = (1/(1+obj.factor*main.params.tau{primalNumber})) * (main.xTilde{primalNumber} + (obj.factor*main.params.tau{primalNumber})*obj.f);
        end
    end
end