%represents base class for terms containing the operator
%K[u,w1,w2] where K=[Dx,-I,0,Dy,0,-I]
%corresponds to three primal variable (u, w1, w2)
classdef basicSecondOrderGradient < basicDualizedOperator
    properties
    end

    methods
        function obj = basicSecondOrderGradient(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            initVar('discretization','forward');

            nPx = prod(dims);
            %K = [Dx,-I,0;Dy,0,-I]
            operatorList{1} = gradientOperator(dims,1,'discretization',discretization);
            operatorList{2} = -identityOperator(nPx);
            operatorList{3} = zeroOperator(nPx);
            operatorList{4} = gradientOperator(dims,2,'discretization',discretization);
            operatorList{5} = zeroOperator(nPx);
            operatorList{6} = -identityOperator(nPx);

            %numPrimals is 3
			obj = obj@basicDualizedOperator(alpha,3,operatorList,varargin);
        end
    end
end
