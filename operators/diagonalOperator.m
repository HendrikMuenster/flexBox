%represents a diagonal matrix
classdef diagonalOperator < basicOperator
    properties
        diagVector
    end
    
    methods
        function obj = diagonalOperator(diagVector,varargin)
            obj.diagVector = diagVector;
        end
        
        function result = mtimes(obj,vector)
            result = obj.diagVector .* vector;
        end
        
        function result = abs(obj)
            result = returnMatrix(obj);
        end
        
        function mat = returnMatrix(obj)
            mat = spdiags(obj.diagVector,0,numel(obj.diagVector),numel(obj.diagVector));
        end
        
%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end
    end
    
end

