%represents base class for terms containing the shear operator
%K(u,v) =
%corresponds to two primal variables (u,v)
classdef basicShear < basicDualizedOperator
    properties
    end

    methods
        function obj = basicShear(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradientND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            end

            for i=1:numel(dims)
                operatorList{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: )';
            end

            operatorList{2} = -operatorList{2};
            operatorList{3} = operatorList{2};
            operatorList{4} = operatorList{1};

            obj = obj@basicDualizedOperator(alpha,numel(operatorList),operatorList,varargin);

        end
    end
end
