%prox for F = alpha / 2 |\nabla u|^2
classdef L1dataTermOperator < basicDualizedDataterm
    methods
        function obj = L1dataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,A,f,varargin);
            
            obj.CPPsupport = 1;
        end
        
        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:)));
        end
        
    end
end