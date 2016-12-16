%represents base class for terms containing the curl operator
%K(u,v) = [u_y^T - v_x^T]
%corresponds to two primal variables (u,v)
classdef basicCurl < basicDualizedOperator
    properties
    end

    methods
        function obj = basicCurl(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradND( dims,ones(numel(dims),1) );
            end

            for i=1:numel(dims)
                operatorList{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
            end

            if (exist('switchedCurl','var') && switchedCurl == 1)
                operatorList{1} = operatorList{1}';
                operatorList{2} = operatorList{2}';
            else
                tmpOp = operatorList{1};
                operatorList{1} = operatorList{2}';
                operatorList{2} = tmpOp';
            end

            if (exist('noMinus','var') && noMinus == 1)
            else
                operatorList{2} = -operatorList{2};
            end

            obj = obj@basicDualizedOperator(alpha,numel(operatorList),operatorList,varargin);
        end
    end
end
