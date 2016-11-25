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
            %vararginParser;

            initVar('stepsize',ones(numel(inputDimension),1));
            initVar('discretization','forward');

            obj.inputDimension = inputDimension;
            obj.gradDirection = gradDirection;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                obj.type = 'backward';
                obj.matrix = generateBackwardGradND( inputDimension,stepsize,gradDirection );
            else
                obj.type = 'forward';
                obj.matrix = generateForwardGradND( inputDimension,stepsize,gradDirection );
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

        function result = size(obj,dim)
            result = size(obj.matrix,dim);
        end

        function result = getMaxRowSumAbs(obj)
            result = 2;
        end
    end

end
