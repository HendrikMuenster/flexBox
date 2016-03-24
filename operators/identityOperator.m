%represents an identity matrix
classdef identityOperator < basicOperator
    properties
        nPx;
    end
    
    methods
        function obj = identityOperator(nPx,varargin)
            obj.nPx = nPx;
        end
        
        function result = mtimes(obj,vector)
            result = vector;
        end
        
        function result = abs(obj)
            result = returnMatrix(obj);
        end
        
        function mat = returnMatrix(obj)
            mat = speye(obj.nPx);
        end
        
%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end
    end
    
end

