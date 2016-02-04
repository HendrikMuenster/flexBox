classdef L2dataTermOperator < basicDualizedDataterm
    methods
        function obj = L2dataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,A,f,varargin);
        end
        
        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = (obj.factor/(main.params.sigma{dualNumber}+obj.factor)) * (main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:));
        end
        
    end
end