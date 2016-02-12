%represents \| K[u,w1,w2] \| with K=[Dx,-I,0,Dy,0,-I]
classdef basicSecondOrderGradient < basicDualizedOperator & tildeMultiOperatorMultiDual
    properties
    end
    
    methods
        function obj = basicSecondOrderGradient(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradientND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            end

            %K = [Dx,-I;Dy,-I]
            operatorList{1} = opTmp( 1 : prod(dims),: );
            operatorList{2} = -speye(prod(dims));
            operatorList{3} = sparse(prod(dims),prod(dims));
            operatorList{4} = opTmp( prod(dims) + 1 : 2 * prod(dims),: );
            operatorList{5} = sparse(prod(dims),prod(dims));
            operatorList{6} = -speye(prod(dims));
            
            %numPrimals is 3
			obj = obj@basicDualizedOperator(alpha,3,operatorList,varargin);
        end
    end
end