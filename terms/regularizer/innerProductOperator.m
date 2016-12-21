%represents the term
%\alpha <b,Au>
%corresponds to one primal variable u and an arbitrary linear operator A
classdef innerProductOperator < basicDualizedDataterm & innerProductProxDual
    methods
        function obj = innerProductOperator(alpha,numPrimals,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,f,varargin);
        end
    end
end
