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
            
            obj.isMinus = 0;
        end

        function result = mtimes(obj,vector)
%             %this is a lightweight test
%             if (strcmp(obj.type,'forward') && obj.gradDirection == 1)
%                 result = imfilter(reshape(vector,obj.inputDimension),[0,-1,1]');
%                 result(end,:) = 0;
%             elseif (strcmp(obj.type,'forward') && obj.gradDirection == 2)
%                 result = imfilter(reshape(vector,obj.inputDimension),[0,-1,1]);
%                 result(:,end) = 0;
%             elseif (strcmp(obj.type,'backward') && obj.gradDirection == 1)
%                 %tic;
%                 vector = reshape(vector,obj.inputDimension);
%                 result = imfilter(vector,[1,-1,0]','same');
%                 result(1,:) = -vector(1,:); 
%                 result(end,:) = vector(end-1,:);% toc;
%                 %tic;result2 = obj.matrix * vector(:);toc;
%                 %result2 = reshape(result2,obj.inputDimension);
%                 
%                 
%                 %figure(1);imagesc(result - result2);
%                 %pause
%             elseif (strcmp(obj.type,'backward') && obj.gradDirection == 2)
%                 vector = reshape(vector,obj.inputDimension);
%                 result = imfilter(vector,[1,-1,0],'same');
%                 result(:,1) = -vector(:,1); 
%                 result(:,end) = vector(:,end-1);
%             else
%                 result = obj.matrix * vector;
%             end
            result = obj.matrix * vector;
            
            if (obj.isMinus)
                result = -result;
            end
            
            result = result(:);
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

        function result = size(obj,varargin)
            if (nargin < 2)
                result = size(obj.matrix);
            else
                result = size(obj.matrix,varargin{1});
            end
        end
        
        function result = getRowSumAbs(obj)
            result = sum(abs(obj.matrix),2);
        end
        
        function result = getMaxRowSumAbs(obj)
            result = max(obj.getRowSumAbs());
        end
    end

end
