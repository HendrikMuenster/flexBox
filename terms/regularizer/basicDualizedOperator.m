%represents the base class for all dualized regularization terms
classdef basicDualizedOperator < basicDualizedDataterm
    methods
        function obj = basicDualizedOperator(alpha,numPrimals,A,varargin)
            if (nargin > 3 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            
            if (iscell(A))
                %A is cell construct create list of fs
                listF = cell(numel(A) / numPrimals,1);
            else
                %A is sparse or full matrix, create zero vector
                listF = {zeros(size(A,1),1)};
            end
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,listF,varargin);

%             obj = obj@dualPart(alpha);
% 
%             if (iscell(A))
%                 obj.numVars = numel(A) / numPrimals;
%                 obj.numPrimals = numPrimals;
% 
% 
%                 for i=1:numPrimals
%                     obj.myTau{i} = 0;
%                 end
%                 for i=1:numel(A) / numPrimals
%                     obj.mySigma{i} = 0;
%                 end
% 
% 
%                 for i=1:numel(A) / numPrimals
%                     for j=1:numPrimals
%                         opNum = (i-1)*numPrimals + j;
% 
%                         opTmp = A{opNum};
%                         if (iscell(opTmp))
%                             opTmp = opTmp{:};%uncell
%                         end
% 
%                         obj.operator{opNum} = opTmp;
%                         obj.operatorT{opNum} = opTmp';
%                         obj.length{opNum} = size(opTmp,1);
% 
%                         obj.mySigma{i} = obj.mySigma{i} + max(sum(abs(opTmp),1));
%                         obj.myTau{j} = obj.myTau{j} + max(sum(abs(opTmp),2));
%                     end
%                 end
%             else
%                 obj.numVars = 1;
%                 obj.length{1} = size(A,1);
%                 obj.operator{1} = A;
%                 obj.operatorT{1} = A';
% 
%                 obj.mySigma{1} = max(sum(abs(A),1));
%                 obj.myTau{1} = max(sum(abs(A),2));
%             end


        end
    end
end
