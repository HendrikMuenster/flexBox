%represents the term
%\alpha |Au|_F (frobenius norm)
%corresponds to one primal variable u and an arbitrary linear operator A
classdef frobeniusOperator < basicDualizedOperator & FrobeniusProxDual
    properties
    end

    methods
        function obj = frobeniusOperator(alpha,numPrimals,A,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
        end
    end
end
