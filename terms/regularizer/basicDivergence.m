%represents base class for terms containing the divergence operator
%\nabla \cdot u
classdef basicDivergence < basicDualizedOperator
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

            opNum = 1;
            for i=1:numel(dims)
                operatorList{opNum} = gradientOperator(dims,i,varargin);
                opNum = opNum + 1;
            end

			%numPrimals is numel(operatorList)
			obj = obj@basicDualizedOperator(alpha,numel(operatorList),operatorList,varargin);
        end
    end
end
