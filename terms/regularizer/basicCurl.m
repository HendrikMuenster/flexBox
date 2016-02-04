%corresponds to two primal variables (u,v) and represents the curl operator operator 
%K(u,v) = [u_y^T - v_x^T]
classdef basicCurl < dualPart & tildeMultiOperatorMultiDual
    properties
    end
    
    methods
        function obj = basicCurl(alpha,dims,varargin)
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
            obj.numVars = 1; %divergence produces scalar quantity
            for i=1:numel(dims)
                obj.length{i} = prod(dims);
                obj.operator{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                obj.myTau{i} = 2;
            end
            obj.mySigma{1} = 2*numel(dims);
            
            if (exist('switchedCurl','var') && switchedCurl == 1)
                obj.operator{1} = obj.operator{1}';
                obj.operator{2} = obj.operator{2}';
            else
                tmpOp = obj.operator{1};
                obj.operator{1} = obj.operator{2}';
                obj.operator{2} = tmpOp';
            end
            
            if (exist('noMinus','var') && noMinus == 1)
            else
                obj.operator{2} = -obj.operator{2};
            end
        end
    end
end