%
classdef basicGradient < basicDualizedOperator & tildeMultiOperatorMultiDual
    methods
        function obj = basicGradient(alpha,dims,varargin)
            if (nargin > 2 && numel(varargin) == 1)
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
            for i=1:numel(usedims)
                if (usedims(i) == 1)
                    operatorList{opNum} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                    opNum = opNum + 1;
                end
            end
			
			%numPrimals is 1
			obj = obj@basicDualizedOperator(alpha,1,operatorList,varargin);
        end


    end
end