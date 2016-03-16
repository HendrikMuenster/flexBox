%class for Kullback-Leibler Divergence, incorporates positivity of Ku
% \min_u Ku-f + f\log(f/Ku) s.t. u>=0
classdef KLdataTermOperator < basicDualizedDataterm
    methods
        function obj = KLdataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,A,f,varargin);
            
            obj.f = max(0,obj.f); %f has to be positive
            
            obj.CPPsupport = 1;
        end
        
        function applyProx(obj,main,dualNumber,~)
            main.y{dualNumber} = 0.5*(1 + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}-1).^2 + 4*main.params.sigma{dualNumber}*obj.f(:) ));
        end
        
    end
end