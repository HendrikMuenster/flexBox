%term for \alpha |Au-f|_1
classdef L1dataTermOperator < basicDualizedDataterm & L1DataProxDual
    methods
        function obj = L1dataTermOperator(alpha,A,f,varargin)
            obj = obj@basicDualizedDataterm(alpha,1,A,f,varargin);
        end
    end
end