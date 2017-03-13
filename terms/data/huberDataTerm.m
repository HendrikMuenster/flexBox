%term for \alpha |u-f|_{H_epsi}
classdef huberDataTerm < huberDataTermOperator & HuberDataProxDual
    methods
        function obj = huberDataTerm(alpha,f,epsi,varargin)
			obj = obj@huberDataTermOperator(alpha,identityOperator(numel(f)),f,epsi);
        end
    end
end