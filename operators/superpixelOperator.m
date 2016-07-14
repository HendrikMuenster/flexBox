%downsamples image of size #factor*#targetDimension to size #targetDimension
%therefore #factor must be a positive integer value
classdef superpixelOperator < basicOperator
    properties
        targetDimension
        factor
        matrix
    end
    
    methods
        
        function obj = superpixelOperator(targetDimension,factor,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            obj.targetDimension = targetDimension;
            obj.factor = factor;
            
            
            %generate matrix
            mask = reshape(1:prod(targetDimension),targetDimension);

            mask = kron(mask,ones(factor));

            weight = 1/factor^2;

            obj.matrix = sparse(mask(:),1:(prod(targetDimension)*factor^2),weight);
        end
        
        function result = mtimes(obj,vector)
            result = obj.matrix * vector;
        end
        
        function result = abs(obj)
            result = abs(obj.matrix);
        end
        
        function mat = returnMatrix(obj)
            mat = obj.matrix;
        end
        
        function res = ctranspose(obj)
            res = obj;
            res.matrix = res.matrix';
        end
        
        function result = size(obj,dim)
            result = size(obj.matrix,dim);
        end
        
        function result = getMaxRowSumAbs(obj)
            result = 1;
        end
    end
    
end

