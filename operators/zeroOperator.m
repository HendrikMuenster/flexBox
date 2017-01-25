%represents an empty matrix
classdef zeroOperator < basicOperator
    properties
        nPx
    end
    
    methods
        function obj = zeroOperator(nPx,varargin)
            obj.nPx = nPx;
        end
        
        function result = mtimes(~,~)
            result = 0;
        end

        function result = abs(obj)
            result = returnMatrix(obj);
        end
        
        function mat = returnMatrix(obj)
            mat = sparse(obj.nPx,obj.nPx);
        end
        
        
        function result = size(obj,varargin)
            result = obj.nPx;
            
            if (nargin < 2)
                result = [result,result]; 
            end
        end
        
%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end

        function result = getMaxRowSumAbs(obj)
            result = 0;
        end
    end
    
end

