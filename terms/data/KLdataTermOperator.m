%class for Kullback-Leibler Divergence, incorporates positivity of Ku
% \min_u Ku-f + f\log(f/Ku) s.t. u>=0
classdef KLdataTermOperator < basicDualizedDataterm & KLDataProxDual
    methods
        function obj = KLdataTermOperator(alpha,A,f,varargin)
            if (iscell(A))
                numPrimals = numel(A);
            else
                numPrimals = 1;
            end
		
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,f,varargin);
            
            obj.f{1} = max(0,obj.f{1}); %f has to be positive
        end
    end
end