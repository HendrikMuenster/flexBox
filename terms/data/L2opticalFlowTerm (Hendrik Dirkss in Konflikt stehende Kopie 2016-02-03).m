%prox for G = alpha / 2 |u_t + \nabla u\cdot v|^2, where v is the unknown
classdef L2opticalFlowTerm < basicOpticalFlow

    methods
        function obj = L2opticalFlowTerm(alpha,image1,image2,varargin)
            obj = obj@basicOpticalFlow(alpha,image1,image2,varargin);
        end
        
        function init(obj,myNumber,main)

        end
        
        function applyProx(obj,main,primalNumbers)
            tau = obj.factor*main.params.tau{primalNumbers(1)};
            
            c1 = 1+tau*obj.ux2;
            c2 = tau*obj.uxuy;
            c3 = 1+tau*obj.uy2;

            teiler = c1.*c3-c2.^2;
            
            b1 = main.xTilde{primalNumbers(1)} - tau*obj.uxut - tau * obj.breg{1};
            b2 = main.xTilde{primalNumbers(2)} - tau*obj.uyut - tau * obj.breg{2};
            
            main.x{primalNumbers(1)} = (b1.*c3-c2.*b2)./teiler;
            main.x{primalNumbers(2)} = (b2.*c1-c2.*b1)./teiler;
        end
    end
end