%represents base class for terms containing the hessian operator
%
%corresponds to one primal variable u
classdef basicHessian < dualPart & tildeMultiOperatorMultiDual
    properties
    end

    methods
        function obj = basicHessian(alpha,dims,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );

            obj = obj@dualPart(alpha);
            obj.numVars = 1;
            obj.length{1} = prod(dims);
            obj.mySigma{1} = 8*numel(dims);
            obj.myTau{1} = 8*numel(dims);

            for i=1:numel(dims)
                if (dims(i) ~= 1)
                    op{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                end
                %obj.operator{i} = op'*op;
            end

            obj.operator{1} = op{1}'*op{1} + op{2}'*op{1} + op{1}'*op{2} + op{2}'*op{2};
            obj.operator{2} = op{1}'*op{2} + op{1}'*op{2} + op{2}'*op{2} + op{1}'*op{1};

        end
    end
end
