%prox for G = alpha  |u-f|
classdef L1dataTerm < L1dataTermOperator
    methods
        function obj = L1dataTerm(alpha,f)
            obj = obj@L1dataTermOperator(alpha,identityOperator(numel(f)),f);
        end
    end
end