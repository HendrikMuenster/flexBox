%represents an identity matrix
classdef identityOperator < basicOperator
    properties
        nPx;
        minus;
    end
    
    methods
        function obj = identityOperator(nPx,varargin)
            obj.nPx = nPx;
            obj.minus = 0;
        end
        
        function result = mtimes(obj,vector)
            if (obj.minus)
                result = -vector;
            else
                result = vector;
            end
        end
        
        function result = abs(obj)
            if (obj.minus)
                result = -returnMatrix(obj);
            else
                result = returnMatrix(obj);
            end
            
        end
        
        function result = returnMatrix(obj)
            if (obj.minus)
                result = -speye(obj.nPx);
            else
                result = speye(obj.nPx);
            end
        end
        
        function result = uminus(obj)
            result = obj;
            result.minus = ~result.minus;
        end
        
        function result = size(obj,dim)
            result = obj.nPx;
        end
        
%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end

        function result = getMaxRowSumAbs(obj)
            result = 1;
        end
    end
    
end

