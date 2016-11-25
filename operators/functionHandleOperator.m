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

        function result = size(obj,dim)
            sizeVec = obj.fHandle(zeros(obj.argumentSize));
            %vectorize result
            sizeVec = size(sizeVec(:));
            result = sizeVec(dim);
        end

        function result = getMaxRowSumAbs(obj)
            result = obj.operatorNorm;
        end
    end

end
