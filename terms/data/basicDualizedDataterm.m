%represents base class for dualized data terms
%all data terms that can be dualized should derive this class directly or indirectly
%a dualized data term is of the form \alpha |Au-f|,
%where A is the linear operator, u is the value to be minimized, and f is the data part
%the class initializes required variables for the fixed-point algorithm
classdef basicDualizedDataterm < dualPart & tildeMultiOperatorMultiDual
    properties
        f;
    end

    methods
        function obj = basicDualizedDataterm(alpha,A,f,varargin)
            if (nargin > 3 && numel(varargin) == 1)
                varargin = varargin{1};
            end

            vararginParser;

            obj = obj@dualPart(alpha);

            %data terms always have only one dual variable
            obj.numVars = 1;

            if (iscell(A))
                obj.numPrimals = numel(A);


                for i=1:obj.numPrimals
                    obj.myTau{i} = 0;
                end
                for i=1:numel(A) / obj.numPrimals
                    obj.mySigma{i} = 0;
                end


                for i=1:numel(A) / obj.numPrimals %for all dual variables
                    for j=1:obj.numPrimals %for all primals
                        opNum = (i-1)*obj.numPrimals + j;

                        opTmp = A{opNum};
                        if (iscell(opTmp))
                            opTmp = opTmp{:};%uncell
                        end

                        obj.operator{opNum} = opTmp;
                        obj.operatorT{opNum} = opTmp';

                        obj.length{opNum} = size(opTmp,1);

                        if (ismatrix(A) || issparse(opTmp))
                            obj.mySigma{i} = obj.mySigma{i} + max(sum(abs(opTmp),1));
                            obj.myTau{j} = obj.myTau{j} + max(sum(abs(opTmp),2));
                        else
                            %this method must be implemented by every
                            %custom operator
                            obj.mySigma{i} = obj.mySigma{i} + obj.operator{opNum}.getMaxRowSumAbs();
                            obj.myTau{j} = obj.myTau{j} + obj.operatorT{opNum}.getMaxRowSumAbs();
                        end
                    end
                end
            else %defaul value numPrimals = 1
                obj.length{1} = size(A,1);
                obj.operator{1} = A;
                obj.operatorT{1} = A';

                if (ismatrix(A) || issparse(A))
                    obj.mySigma{1} = max(sum(abs(obj.operator{1}),1));
                    obj.myTau{1} = max(sum(abs(obj.operator{1}),2));
                else
                    %this method must be implemented by every
                    %custom operator
                    obj.mySigma{1} = obj.operator{1}.getMaxRowSumAbs();
                    obj.myTau{1} = obj.operatorT{1}.getMaxRowSumAbs();
                end
            end

            numRowsOperator = size(obj.operator{1},1);

            obj.f = f(:);

            if (numel(f) ~= numRowsOperator)
                error('Input data f does not fit the number of rows in the operator(s)');
            end

        end
    end
end
