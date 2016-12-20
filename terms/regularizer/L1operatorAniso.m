%represents the term
%\alpha |Au|_1,1
%in the L_1,1 norm
%corresponds to one primal variable u and an arbitrary linear operator A
classdef L1operatorAniso < basicDualizedOperator & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1operatorAniso(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
