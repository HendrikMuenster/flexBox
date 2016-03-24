%represents an empty matrix
classdef zeroOperator < basicOperator
    properties
        nPx
    end
    
    methods
        function obj = zeroOperator(nPx,varargin)
            obj.nPx = nPx;
        end
        
        function result = mtimes(obj,vector)
            result = 0;
        end

        function result = abs(obj)
            result = returnMatrix(obj);
        end
        
        function mat = returnMatrix(obj)
            mat = sparse(obj.nPx,obj.nPx);
        end
        
%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end
    end
    
end

