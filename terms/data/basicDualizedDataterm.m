%
classdef basicDualizedDataterm < dualPart & tildeMultiOperatorMultiDual
    properties
        f;
    end
    
    methods
        function obj = basicDualizedDataterm(alpha,A,f,varargin)
            if (nargin > 3 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            obj = obj@dualPart(alpha);
            obj.numVars = 1;
            obj.length{1} = size(A,1);
            obj.operator{1} = A;
            obj.f = f(:);
            obj.mySigma{1} = max(sum(abs(A),1));
            obj.myTau{1} = max(sum(abs(A),2));
        end
    end
end