%represents the concatenation of two operators A and B
classdef concatOperator < basicOperator
    properties
        A;          %operator for A
        B;          %operator for B
        operation   %operation to represent (composition, addition, difference)
    end

    methods
        function obj = concatOperator(A,B,operation,varargin)
            obj.A = A;
            obj.B = B;
            
            if ((strcmp(operation,'addition') || strcmp(operation,'difference')) && size(A,1) ~= size(B,1))
                error('Cannot create concat operator, because size(A,1) ~= size(B,1)');
            end
            
            
            if (strcmp(operation,'composition') && size(A,2) ~= size(B,1))
                error('Cannot create concat operator, because size(A,2) ~= size(B,1)');
            end
            
            if (strcmp(operation,'composition') || strcmp(operation,'addition') || strcmp(operation,'difference'))
                obj.operation = operation;
            else
                error('Unknown operation for concat operator. Please choose between "composition", "addition" and "difference".');
            end
        end

        function result = mtimes(obj,vector)
            if (strcmp(obj.operation,'composition'))
                result = obj.A * (obj.B * vector(:));
            elseif (strcmp(obj.operation,'addition'))
                result = obj.A * vector(:) + obj.B * vector(:);
            elseif (strcmp(obj.operation,'difference'))
                result = obj.A * vector(:) - obj.B * vector(:);
            end
            
            if (obj.isMinus)
                result = -result;
            end
        end

        function result = size(obj,varargin)
            if (nargin < 2)
                if (strcmp(obj.operation,'composition'))
                    result = [size(obj.A,1),size(obj.B,1)];
                elseif (strcmp(obj.operation,'addition') || strcmp(obj.operation,'difference'))
                    result = [size(obj.A,1),size(obj.A,2)];
                end
            else
                dim = varargin{1};
                if (dim == 1)
                    result = size(obj.A,1);
                else
                    result = size(obj.A,2);
                end
            end
        end

        function result = returnMatrix(obj)
            if (strcmp(obj.operation,'composition'))
                result = obj.A.returnMatrix() * obj.B.returnMatrix();
            elseif (strcmp(obj.operation,'addition') || strcmp(obj.operation,'difference'))
                result = obj.A.returnMatrix() + obj.B.returnMatrix();
            end
        end

        function res = ctranspose(obj)
            res = obj;
            if (strcmp(obj.operation,'composition'))
                tmpOp = res.A;
                res.A = res.B';
                res.B = tmpOp';
            elseif (strcmp(obj.operation,'addition') || strcmp(obj.operation,'difference'))
                res.A = res.A';
                res.B = res.B';
            end
        end
        
        function result = abs(obj)
            if (strcmp(obj.operation,'composition'))
                result = abs(obj.A) * abs(obj.B);
            elseif (strcmp(obj.operation,'addition') || strcmp(obj.operation,'difference'))
                result = abs(obj.A) + abs(obj.B);
            end
        end

        function result = getRowSumAbs(obj)
            if (strcmp(obj.operation,'composition'))
                if (isa(obj.A,'basicOperator'))
                    resultA = obj.A.getRowSumAbs();
                else
                    resultA = (sum(abs(obj.A),2));
                end
                if (isa(obj.B,'basicOperator'))
                    resultB = obj.B.getRowSumAbs();
                else
                    resultB = (sum(abs(obj.B),1));
                end
                
                result = max(resultA(:)) * max(resultB(:));
            elseif (strcmp(obj.operation,'addition') || strcmp(obj.operation,'difference'))
                if (isa(obj.A,'basicOperator'))
                    resultA = obj.A.getRowSumAbs();
                else
                    resultA = (sum(abs(obj.A),2));
                end
                if (isa(obj.B,'basicOperator'))
                    resultB = obj.B.getRowSumAbs();
                else
                    resultB = (sum(abs(obj.B),2));
                end
                
                result = resultA + resultB;
            end
        end
        
        function result = getMaxRowSumAbs(obj)
            result = getRowSumAbs(obj);
        end
        
    end

end
