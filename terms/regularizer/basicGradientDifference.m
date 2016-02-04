%
classdef basicGradientDifference < dualPart & tildeMultiOperatorMultiDual
    methods
        function obj = basicGradientDifference(alpha,dims,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            obj = obj@dualPart(alpha);
            obj.numPrimals = 2; %expecting two primal variables
            
            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradientND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            end
            
            %usedims should be a {0,1} array of length dims indicating whether a
            %dimension should be used or not
            if (exist('usedims','var'))
                obj.numVars = sum(usedims);
            else
                obj.numVars = numel(dims);
                usedims = ones(obj.numVars,1);
            end
            
            opNum = 1;
            for i=1:numel(usedims)
                if (usedims(i) == 1)
                    obj.length{opNum} = prod(dims);
                    obj.length{opNum+1} = prod(dims);
                    obj.operator{opNum} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: );
                    obj.operator{opNum+1} = -obj.operator{opNum};
                    obj.mySigma{opNum} = 4;
                    obj.mySigma{opNum+1} = 4;
                    opNum = opNum + 2;
                end
            end
            obj.myTau{1} = 2*sum(usedims);
            obj.myTau{2} = 2*sum(usedims);
        end


    end
end