%%represents base class for terms containing the laplace operator
%\Delta u
%corresponds to one primal variable u
classdef basicLaplace < basicDualizedOperator
    properties
    end

    methods
        function obj = basicLaplace(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );

            for i=1:numel(dims)
                op = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                operatorListTmp{i} = op'*op;
            end

            operatorList{1} = sparse(prod(dims),prod(dims));
            for i=2:numel(dims)
                operatorList{1} = operatorList{1} + operatorListTmp{i};
            end

            obj = obj@basicDualizedOperator(alpha,1,operatorList,varargin);
        end
    end
end
