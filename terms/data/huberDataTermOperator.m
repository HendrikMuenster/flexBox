%term for \alpha |Au-f|_{H_epsi}
classdef huberDataTermOperator < basicDualizedDataterm & HuberDataProxDual
    properties
        epsi
    end
    
    methods
        function obj = huberDataTermOperator(alpha,A,f,epsi,varargin)
            if (iscell(A))
                numPrimals = numel(A);
            else
                numPrimals = 1;
            end
		
            obj = obj@basicDualizedDataterm(alpha,numPrimals,A,f,varargin);
            obj.epsi = epsi;
        end
    end
end