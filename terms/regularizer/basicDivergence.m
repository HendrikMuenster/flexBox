%
classdef basicDivergence < basicDualizedOperator & tildeMultiOperatorMultiDual
    properties
    end
    
    methods
        function obj = basicDivergence(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            %usedims should be a {0,1} array of length dims indicating whether a
            %dimension should be used or not
            if (~exist('usedims','var'))
                usedims = ones(numel(dims),1);
            end
            
            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradientND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            end
            
            opNum = 1;
            for i=1:numel(dims)
                if (usedims(i) == 1 && dims(i) ~= 1)
                    operatorList{opNum} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                    opNum = opNum + 1;
                end
            end

			%numPrimals is numel(operatorList)
			obj = obj@basicDualizedOperator(alpha,numel(operatorList),operatorList,varargin);
        end
    end
end