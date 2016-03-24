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

            nPx = prod(dims);
            %K = [Dx,-I,0;Dy,0,-I]
            operatorList{1} = opTmp( 1 : nPx,: );
            operatorList{2} = -identityOperator(nPx);
            operatorList{3} = zeroOperator(nPx);
            operatorList{4} = opTmp( nPx + 1 : 2 * nPx,: );
            operatorList{5} = zeroOperator(nPx);
            operatorList{6} = -identityOperator(nPx);
            
            %numPrimals is 3
			obj = obj@basicDualizedOperator(alpha,3,operatorList,varargin);
        end
    end
end