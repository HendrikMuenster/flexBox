classdef basicIdentity < dualPart & tildeSingleOperator
    properties
    end
    
    methods
        function obj = basicIdentity(alpha,dims,varargin)
            obj = obj@dualPart(alpha);
            obj.numVars = 1;
            obj.length{1} = prod(dims);
            obj.mySigma{1} = 1;
            obj.myTau{1} = 1;
            obj.operator{1} = speye(prod(dims));
        end
    end
end