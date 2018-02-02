% represents a subsampling Operator
% acts like permutationOperator but allows repetition of elements
% expects a vector of indices (ind) and effectivley performs 
% vec(ind) in Operator Vector Multiplication
classdef subsamplingOperator < basicOperator
    properties
        indices  %vector/matrix of indices
        numElem  %number of elements in input vector
        matrix   %matrix representation of operator
    end
    
    methods (Access = private)
        function mat = generateMatrix(obj)
            rows = (1:length(obj.indices))';
            mat = sparse(rows, obj.indices, true(length(obj.indices),1), length(obj.indices), obj.numElem);
        end
    end

    methods
        function obj = subsamplingOperator(indices, inputDimension, varargin)
            obj.indices = indices(:);
            obj.numElem = prod(inputDimension);
            obj.isMinus = 0;
            obj.matrix = obj.generateMatrix();
        end

        function result = mtimes(obj,vector)
            result = obj.matrix * vector(:);

            if (obj.isMinus)
                result = -result;
            end
            
            result = result(:);
        end

        function result = abs(obj)
            result = returnMatrix(obj);
        end

        function mat = returnMatrix(obj)
            mat = obj.matrix;
        end

        function result = size(obj,varargin)
            opSize = [length(obj.indices), obj.numElem];
            if (nargin > 1)
                dim = varargin{1};
                result = opSize(dim);
            else
                result = opSize;
            end
        end

         function res = ctranspose(obj)
             res = obj;
             res.matrix = res.matrix';
         end
         
        function result = getRowSumAbs(obj)
            result = sum(abs(obj.matrix),2);
        end
        
        function result = getMaxRowSumAbs(obj)
            result = max(obj.getRowSumAbs());
        end
    end

end
