%represents the term
%\alpha |Au|_1,2
%in the L_1,2 norm
%corresponds to one primal variable u and an arbitrary linear operator A
classdef L1operatorIso < basicDualizedOperator & L1IsoProxDual
    properties
    end

    methods
        function obj = L1operatorIso(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
