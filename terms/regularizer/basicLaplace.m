classdef basicLaplace < dualPart & tildeMultiOperatorMultiDual
    properties
    end
    
    methods
        function obj = basicLaplace(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            
            obj = obj@dualPart(alpha);
            obj.numVars = 1; %divergence produces scalar quantity
            obj.length{1} = prod(dims);
            obj.mySigma{1} = 3*numel(dims);
            obj.myTau{1} = 4*numel(dims);
            
            for i=1:numel(dims)
                op = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                obj.operator{i} = op'*op;
            end
            
            for i=2:numel(dims)
                obj.operator{1} = obj.operator{1} + obj.operator{i};
                obj.operator{i} = []; %save memory
            end
        end
    end
end