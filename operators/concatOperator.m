%represents the concatenation of two operators A and B
classdef concatOperator < basicOperator
    properties
        A;          %operator for A
        B;          %operator for B
        transposed; %is-transposed-flag in order to decide which operators to use
    end

    methods
        function obj = concatOperator(A,B,varargin)
            obj.A = A;
            obj.B = B;
            obj.transposed = 0;
        end

        function result = mtimes(obj,vector)
            result = obj.A * (obj.B * vector(:));
            
            if (obj.isMinus)
                result = -result;
            end
        end

        function result = abs(obj)
            result = abs(obj.A) * abs(obj.B);
        end

        function result = size(obj,varargin)
            if (nargin < 2)
                result = [size(obj.A,1),size(obj.B,1)];
            else
                dim = varargin{1};
                if (dim == 1)
                    result = size(obj.A,1);
                else
                    result = size(obj.B,1);
                end
            end
        end

        function result = returnMatrix(obj)
            result = obj.A.returnMatrix() * obj.B.returnMatrix();
        end

        function res = ctranspose(obj)
            res = obj;
            
            tmpOp = res.A;
            res.A = res.B';
            res.B = tmpOp';
        end

        function result = getMaxRowSumAbs(obj)
            if (issparse(obj.A))
                resultA = max(sum(abs(obj.A),1));
            else
                resultA = obj.A.getMaxRowSumAbs();
            end
            if (issparse(obj.B))
                resultB = max(sum(abs(obj.B),1));
            else
                resultB = obj.B.getMaxRowSumAbs();
            end

            result = resultA * resultB;
        end
    end

end
