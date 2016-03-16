classdef L1dataTermOperator < basicDualizedDataterm & L1DataProxDual
    methods
        function obj = L1dataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,A,f,varargin);
            
            obj.CPPsupport = 1;
        end
    end
end