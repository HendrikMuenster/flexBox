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
            
            checkF = sum(obj.f{1}<=0);
            if (checkF > 0)
                disp('For KLdataTerm f has to be strictly positive! Setting critical values to 1e-6');
                obj.f{1} = max(1e-6,obj.f{1});
            end
        end
    end
end