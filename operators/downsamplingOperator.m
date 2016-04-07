%downsamples image of size #inputDimension to size #outputDimension
classdef downsamplingOperator < basicOperator
    properties
        inputDimension
        outputDimension
        matrix
    end
    
    methods
        
        function obj = downsamplingOperator(inputDimension,outputDimension,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            obj.inputDimension = inputDimension;
            obj.outputDimension = outputDimension;
            obj.matrix = generateDownsamplingMatrix( outputDimension,inputDimension );
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

