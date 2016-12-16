%term for \alpha/2 |Au-f|_2^2
classdef L2dataTermOperator < basicDualizedDataterm & L2DataProxDual
    methods
        function obj = L2dataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,1,A,f,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end