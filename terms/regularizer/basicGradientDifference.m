%represents base class for terms containing the gradient differences
%\nabla(u-w)
%corresponds to two primal variables (u,v)
%TODO: generalize for higher dimensions
classdef basicGradientDifference < basicDualizedOperator
    methods
        function obj = basicGradientDifference(alpha,dims,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradND( dims,ones(numel(dims),1) );
            end

            %usedims should be a {0,1} array of length dims indicating whether a
            %dimension should be used or not
            if (~exist('usedims','var'))
                usedims = ones(2,1);
            end

            opNum = 1;
            for i=1:numel(usedims)
                if (usedims(i) == 1 && dims(i) ~= 1)
                    operatorList{opNum} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                    operatorList{opNum+1} = -operatorList{opNum};

                    opNum = opNum + 2;
                end
            end

            obj = obj@basicDualizedOperator(alpha,2,operatorList,varargin);
        end
    end
end
