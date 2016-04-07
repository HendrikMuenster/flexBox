%represents the concatenation of two operators A and B
classdef concatOperator < basicOperator
    properties
        A;
        B;
        AT;
        BT;
        transposed;
    end
    
    methods
        function obj = concatOperator(A,B,varargin)
            obj.A = A;
            obj.B = B;
            obj.AT = A';
            obj.BT = B';
            obj.transposed = 0;
        end
        
        function result = mtimes(obj,vector)
            if (obj.transposed)
                result = obj.BT * (obj.AT * vector);
            else
                result = obj.A * (obj.B * vector);
            end
        end
        
        function result = abs(obj)
            if (obj.transposed)
                result = abs(obj.BT) * abs(obj.AT);
            else
                result = abs(obj.A) * abs(obj.B);
            end
        end
        
        function result = size(obj,varargin)
            
            if (nargin < 2)
                if (obj.transposed)
                    result = [size(obj.B,2),size(obj.A,2)];
                else
                    result = [size(obj.A,1),size(obj.B,1)];
                end
                
            else
                dim = varargin{1};
                if (dim == 1)
                    if (obj.transposed)
                        result = size(obj.B,2);
                    else
                        result = size(obj.A,1);
                    end
                else
                    if (obj.transposed)
                        result = size(obj.A,2);
                    else
                        result = size(obj.B,1);
                    end
                end
            end
        end
        
        function result = returnMatrix(obj)
            result = obj.A.returnMatrix() * obj.B.returnMatrix();
        end
        
        function res = ctranspose(obj)
            res = obj;
            if (res.transposed == 1)
                res.transposed = 0;
            else
                res.transposed = 1;
            end
        end
        
        function result = getMaxRowSumAbs(obj)
            if (obj.transposed == 1)
                result = obj.BT.getMaxRowSumAbs() * obj.AT.getMaxRowSumAbs();
            else
                result = obj.A.getMaxRowSumAbs() * obj.B.getMaxRowSumAbs();
            end
        end
    end
    
end

