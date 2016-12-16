%class for Kullback-Leibler Divergence, incorporates positivity of Ku
% \min_u Ku-f + f\log(f/Ku) s.t. u>=0
classdef KLdataTermOperator < basicDualizedDataterm & KLDataProxDual
    methods
        function obj = KLdataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,1,A,f,varargin);
            
            obj.f{1} = max(0,obj.f{1}); %f has to be positive
            
            obj.CPPsupport = 1;
        end
    end
end