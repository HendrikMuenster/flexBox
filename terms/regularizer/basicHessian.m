%represents base class for terms containing the hessian operator
%
%corresponds to one primal variable u
classdef basicHessian < basicDualizedOperator
    properties
    end

    methods
        function obj = basicHessian(alpha,dims,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );

            for i=1:numel(dims)
                if (dims(i) ~= 1)
                    op{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                end
            end

            operatorList{1} = op{1}'*op{1} + op{2}'*op{1} + op{1}'*op{2} + op{2}'*op{2};
            operatorList{2} = op{1}'*op{2} + op{1}'*op{2} + op{2}'*op{2} + op{1}'*op{1};

            obj = obj@basicDualizedOperator(alpha,numel(operatorList),operatorList,varargin);
        end
    end
end
