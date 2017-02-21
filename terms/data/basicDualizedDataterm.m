%represents base class for dualized data terms
%all data terms that can be dualized should derive this class directly or indirectly
%a dualized data term is of the form \alpha |Au-f|,
%where A is the linear operator, u is the value to be minimized, and f is the data part
%the class initializes required variables for the fixed-point algorithm
classdef basicDualizedDataterm < dualPart
    properties
        f;%cell array of data
    end

    methods
        function obj = basicDualizedDataterm(alpha,numPrimals,A,f,varargin)
            if (nargin > 4 && numel(varargin) == 1)
                varargin = varargin{1};
            end

            vararginParser;

            obj = obj@dualPart(alpha);
            
            if (~iscell(A))
                A = {A};
            end
            
            if (~iscell(f))
                f = {f};
            end

            obj.numVars = numel(A) / numPrimals;
            obj.numPrimals = numPrimals;

            for i=1:obj.numVars %for all dual variables
                if (numel(f{i}) > 0)
                    obj.f{i} = reshape(f{i},[numel(f{i}),1]);
                else
                    obj.f{i} = {};
                end
                
                for j=1:obj.numPrimals %for all primals
                    opNum = (i-1)*obj.numPrimals + j;

                    opTmp = A{opNum};
                    if (iscell(opTmp))
                        opTmp = opTmp{:};%uncell
                    end

                    if (numel(obj.f{i}) ~= size(opTmp,1) && numel(obj.f{i}) > 0)
                        error(['Input data f{',num2str(i),'} does not fit the number of rows in the operator ',opNum]);
                    end

                    obj.operator{opNum} = opTmp;
                    obj.operatorT{opNum} = opTmp';

                    obj.length{opNum} = size(opTmp,1);
                end
            end
        end
        
        function yTilde(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                main.yTilde{dualNumbers(i)} = main.y{dualNumbers(i)};
            end

            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.yTilde{dualNumbers(i)} = main.yTilde{dualNumbers(i)} + main.params.sigma{dualNumbers(i)} * (obj.operator{operatorNumber}* main.xBar{primalNumbers(j)});
                end
            end
        end

        function xTilde(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.xTilde{primalNumbers(j)} = main.xTilde{primalNumbers(j)} - main.params.tau{primalNumbers(j)}*(obj.operatorT{operatorNumber} * main.y{dualNumbers(i)});
                end
            end
        end

        function yError(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.yError{dualNumbers(i)} = main.yError{dualNumbers(i)} - obj.operator{operatorNumber}* (main.x{primalNumbers(j)}-main.xOld{primalNumbers(j)});
                end
            end
        end

        function xError(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.xError{primalNumbers(j)} = main.xError{primalNumbers(j)} - obj.operatorT{operatorNumber} * (main.y{dualNumbers(i)} - main.yOld{dualNumbers(i)});
                end
            end
        end
    end
end
