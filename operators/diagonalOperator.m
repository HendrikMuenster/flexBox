%represents a diagonal matrix
classdef diagonalOperator < basicOperator
    properties
        diagonalElements  %vector of diagonal elements
    end

    methods
        function obj = diagonalOperator(diagonalElements,varargin)
            obj.diagonalElements = diagonalElements(:);
        end

        function result = mtimes(obj,vector)
            result = obj.diagonalElements .* vector;
            
            if (obj.isMinus)
                result = -result;
            end
        end

        function result = abs(obj)
            result = abs(returnMatrix(obj));
        end

        function mat = returnMatrix(obj)
            mat = spdiags(obj.diagonalElements,0,numel(obj.diagonalElements),numel(obj.diagonalElements));
        end

        function result = size(obj,varargin)
            result = numel(obj.diagonalElements);

            if (nargin < 2)
                result = [result,result];
            end
        end

%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end

        function result = getMaxRowSumAbs(obj)
            result = max(abs(obj.diagonalElements));
        end
    end

end
