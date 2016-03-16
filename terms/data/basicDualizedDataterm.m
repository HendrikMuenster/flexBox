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
%             vararginParser;
%             
%             obj = obj@dualPart(alpha);
%             obj.numVars = 1;
%             obj.length{1} = size(A,1);
%             obj.operator{1} = A;
%             
%             obj.mySigma{1} = max(sum(abs(A),1));
%             obj.myTau{1} = max(sum(abs(A),2));
%             
%             obj.f = f(:);
            
            vararginParser;
            
            obj = obj@dualPart(alpha);
            
            obj.numVars = 1;
            
            if (iscell(A))
                obj.numPrimals = numel(A);
                
               
                for i=1:obj.numPrimals
                    obj.myTau{i} = 0;
                end
                for i=1:numel(A) / obj.numPrimals
                    obj.mySigma{i} = 0;
                end
                    
                
                for i=1:numel(A) / obj.numPrimals
                    for j=1:obj.numPrimals
                        opNum = (i-1)*obj.numPrimals + j;
                        
                        opTmp = A{opNum};
                        if (iscell(opTmp))
                            opTmp = opTmp{:};%uncell
                        end
                        
                        obj.operator{opNum} = opTmp;
                        obj.length{opNum} = size(opTmp,1);
                        
                        obj.mySigma{i} = obj.mySigma{i} + max(sum(abs(opTmp),1));
                        obj.myTau{j} = obj.myTau{j} + max(sum(abs(opTmp),2));
                    end
                end
            else
                obj.length{1} = size(A,1);
                obj.operator{1} = A;

                obj.mySigma{1} = max(sum(abs(A),1));
                obj.myTau{1} = max(sum(abs(A),2));
            end
            
            obj.f = f(:);
            
        end
    end
end