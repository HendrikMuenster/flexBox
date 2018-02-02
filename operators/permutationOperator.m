% represents a permutation Operator
% expects a vector of indices (ind) and effectivley performs 
% vec(ind) in Operator Vector Multiplication
% ind has to contain every 
classdef permutationOperator < basicOperator
    properties
        indices  %vector/matrix of indices
        numElem  %number of elements in input vector
    end
    
    methods (Access = private)
        function mat = generateMatrix(obj)
            rows = (1:length(obj.indices))';
            mat = sparse(rows, obj.indices, true(obj.numElem,1), obj.numElem, obj.numElem);
        end
    end

    methods
        function obj = permutationOperator(indices, varargin)
            obj.indices = indices(:);
            obj.numElem = length(indices);
            obj.isMinus = 0;
            
            if length(indices) ~= length(unique(indices))
                error('Index in indicies was repeated or omitted. Consider using subsamplingOperator');
            end
        end

        function result = mtimes(obj,vector)
            result = vector(obj.indices);

            if (obj.isMinus)
                result = -result;
            end
            result = result(:);
        end

        function result = abs(obj)
            result = obj.generateMatrix();
        end

        function mat = returnMatrix(obj)
            mat = obj.generateMatrix();
        end

        function result = size(obj,varargin)
            opSize = [obj.numElem, obj.numElem];
            if (nargin > 1)
                dim = varargin{1};
                result = opSize(dim);
            else
                result = opSize;
            end
        end

         function res = ctranspose(obj)
             res = obj;
             indicesTmp = [obj.indices (1:length(obj.indices))'];
             indicesTmp = sortrows(indicesTmp);
             res.indices = indicesTmp(:,2);
         end
         
        function result = getRowSumAbs(obj)
            result = ones(obj.numElem, 1);
        end
        
        function result = getMaxRowSumAbs(obj)
            result = 1;
        end
    end
    

end
