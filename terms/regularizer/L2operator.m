%represents the term
%\alpha/2 |Au|_2^2
%in the L_2^2norm
%corresponds to one primal variable u and an arbitrary linear operator A
classdef L2operator < basicDualizedOperator & L2proxDual
    properties
    end

    methods
        function obj = L2operator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
