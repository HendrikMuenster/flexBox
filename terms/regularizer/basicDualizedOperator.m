%
classdef basicDualizedOperator < dualPart & tildeMultiOperatorMultiDual
    methods
        function obj = basicDualizedOperator(alpha,numPrimals,A,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            obj = obj@dualPart(alpha);
            
            if (iscell(A))
                obj.numVars = numel(A) / numPrimals;
                obj.numPrimals = numPrimals;
                
               
                for i=1:obj.numVars
                    obj.myTau{i} = 0;
                end
                for i=1:numel(A) / numPrimals
                    obj.mySigma{i} = 0;
                end
                    
                
                for i=1:numel(A) / numPrimals
                    for j=1:numPrimals
                        opNum = (i-1)*numPrimals + j;
                        
                        obj.operator{opNum} = A{opNum};
                        obj.length{opNum} = size(A{opNum},1);
                        
                        obj.mySigma{i} = obj.mySigma{i} + max(sum(abs(A{opNum}),1));
                        obj.myTau{j} = obj.myTau{i} + max(sum(abs(A{opNum}),2));
                    end
                end
            else
                obj.numVars = 1;
                obj.length{1} = size(A,1);
                obj.operator{1} = A;

                obj.mySigma{1} = max(sum(abs(A),1));
                obj.myTau{1} = max(sum(abs(A),2));
            end
            

        end
    end
end