%term representing alpha |u_t + \nabla u\cdot v|_1, where v is the unknown
classdef L1opticalFlowTerm < basicOpticalFlow & L1DataProxDual

    methods
        function obj = L1opticalFlowTerm(alpha,image1,image2,varargin)
            obj = obj@basicOpticalFlow(alpha,image1,image2,varargin);
            
            obj.CPPsupport = 1;
        end
        
%         function applyProx(obj,main,primalNumbers)
%             tau = obj.factor*main.params.tau{primalNumbers(1)};
%             
%             rho = obj.ut + obj.ux.*main.xTilde{primalNumbers(1)} + obj.uy.*main.xTilde{primalNumbers(2)};
%             
%             c1 = rho<= -tau * obj.normSquared;
%             c2 = rho>= tau * obj.normSquared;
%             
%             main.x{primalNumbers(1)} = main.xTilde{primalNumbers(1)} + tau*obj.ux .* (c1-c2) - obj.nablaNorm{1} .* rho .* (1-c1-c2);
%             main.x{primalNumbers(2)} = main.xTilde{primalNumbers(2)} + tau*obj.uy .* (c1-c2) - obj.nablaNorm{2} .* rho .* (1-c1-c2);
%         end
    end
end