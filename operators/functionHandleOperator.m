% support class for function handles as operators
% needs the operator and adjoint operator as function handle, the size of
% the input argument and the operator norm
classdef functionHandleOperator < basicOperator
    properties
        fHandle
        fHandleT
        argumentSize
        transposed
        operatorNorm
    end

    methods

        function obj = functionHandleOperator(inputfHandle,inputfHandleT,inputArgumentSize,inputOperatorNorm,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            obj.fHandle = inputfHandle;
            obj.fHandleT = inputfHandleT;
            obj.transposed = 0;
            if (isscalar(inputArgumentSize))
                inputArgumentSize = [inputArgumentSize,1];
            end
            
            obj.argumentSize = inputArgumentSize;
            obj.operatorNorm = inputOperatorNorm;

        end

        function result = mtimes(obj,vector)
            if (obj.transposed)
                result = obj.fHandleT(vector);
            else
                result = obj.fHandle(vector);
            end
            result = result(:);
            
            if (obj.isMinus)
                result = -result;
            end
        end

        %this is not correct!
        function value = abs(obj)
            value = obj.operatorNorm;
        end

        function mat = returnMatrix(obj)
            error('returnMatrix not supported');
        end

        function res = ctranspose(obj)
            res = obj;
            res.transposed = ~res.transposed;
        end

        function result = size(obj,varargin)
            sizeVec = obj.fHandle(zeros(obj.argumentSize));
            %vectorize result
            sizeVec = size(sizeVec(:));
            
            if (nargin > 1)
                result = sizeVec(varargin{1});
            else
                result = sizeVec;
            end
        end

        function result = getRowSumAbs(obj)
            result = obj.operatorNorm;
        end
        
        function result = getMaxRowSumAbs(obj)
            result = max(obj.getRowSumAbs());
        end
    end

end
