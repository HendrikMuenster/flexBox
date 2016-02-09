%represents \| K[u,w1,w2] \| with K=[Dx,-I,0,Dy,0,-I]
classdef basicSecondOrderGradient < dualPart & tildeMultiOperatorMultiDual
    properties
    end
    
    methods
        function obj = basicSecondOrderGradient(alpha,dims,varargin)
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
            obj.numPrimals = 3;
            obj.numVars = numel(dims);
            
            %K = [Dx,-I;Dy,-I]
            obj.operator{1} = opTmp( 1 : prod(dims),: );
            obj.operator{2} = -speye(prod(dims));
            obj.operator{3} = sparse(prod(dims),prod(dims));
            obj.operator{4} = opTmp( prod(dims) + 1 : 2 * prod(dims),: );
            obj.operator{5} = sparse(prod(dims),prod(dims));
            obj.operator{6} = -speye(prod(dims));
            
            obj.mySigma{1} = 3;
            obj.mySigma{2} = 3;
            
            obj.myTau{1} = 4;
            obj.myTau{2} = 1;
            obj.myTau{3} = 1;
            
            obj.length{1} = prod(dims);
            obj.length{2} = prod(dims);
        end
    end
end