%represents base class for terms containing the operator
%K(u,v) = [u_x + v_y;u_y - v_x]
%corresponds to two primal variables (u,v) and
classdef basicDivCurl < basicDualizedOperator
    properties
    end

    methods
        function obj = basicDivCurl(alpha,dims,varargin)
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

            operatorList{3} = operatorList{2};
            operatorList{4} = -operatorList{1};

            obj = obj@basicDualizedOperator(alpha,numel(operatorList),operatorList,varargin);

        end


    end
end
