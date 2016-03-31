%represents a gradient operator
classdef gradientOperator < basicOperator
    properties
        matrix
        inputDimension
        gradDirection
        type
    end
    
    methods
        function obj = gradientOperator(inputDimension,gradDirection,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            obj.inputDimension = inputDimension;
            obj.gradDirection = gradDirection;
            
            if (exist('discretization','var') && strcmp(discretization,'backward'))
                obj.type = 'backward';
                obj.matrix = generateBackwardGradientND( inputDimension,ones(numel(inputDimension),1),gradDirection );
            else
                obj.type = 'forward';
                obj.matrix = generateForwardGradientND( inputDimension,ones(numel(inputDimension),1),gradDirection );
            end
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
            
            if (strcmp(res.type,'forward'))
                res.type = 'backward';
            elseif (strcmp(res.type,'backward'))
                res.type = 'forward';
            end
        end
    end
    
end

