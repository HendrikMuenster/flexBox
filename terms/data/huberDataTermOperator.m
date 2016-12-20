%term for \alpha |Au-f|_H_eps
classdef huberDataTermOperator < basicDualizedDataterm & L1HuberDataProxDual
    properties
        epsi
    end
    
    methods
        function obj = huberDataTermOperator(alpha,A,f,epsi,varargin)
            obj = obj@basicDualizedDataterm(alpha,1,A,f,varargin);
            obj.epsi = epsi;
        end
    end
end