%represents base class for terms containing the gradient operator
%\nabla u
%corresponds to one primal variable u
classdef basicGradient < basicDualizedOperator
    methods
        function obj = basicGradient(alpha,dims,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end

            vararginParser;

            %usedims should be a {0,1} array of length dims indicating whether a
            %dimension should be used or not
            initVar('usedims',ones(numel(dims),1));

            opNum = 1;
            for i=1:numel(usedims)
                %if dims(i) equals 1 then matrix is empty
                if (usedims(i) == 1 && dims(i) ~= 1)
                    operatorList{opNum} = gradientOperator(dims,i,varargin);
                    opNum = opNum + 1;
                end
            end

			%numPrimals is 1, numer of duals equals number of gradient directions
			obj = obj@basicDualizedOperator(alpha,1,operatorList,varargin);
        end


    end
end
