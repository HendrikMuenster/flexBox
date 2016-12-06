%represents base class for terms containing the shear operator
%K(u,v) = 
%corresponds to two primal variables (u,v)
classdef basicShear < dualPart & tildeMultiOperatorMultiDual
    properties
    end

    methods
        function obj = basicShear(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradientND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            end

            obj = obj@dualPart(alpha);
            obj.numVars = 2; %divergence produces scalar quantity
            for i=1:numel(dims)
                obj.length{i} = prod(dims);
                obj.operator{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: )';
                obj.myTau{i} = 2;
            end
            obj.mySigma{1} = 3*numel(dims);
            obj.mySigma{2} = 3*numel(dims);

            obj.operator{2} = -obj.operator{2};

            obj.operator{3} = obj.operator{2};
            obj.operator{4} = obj.operator{1};
        end
    end
end
