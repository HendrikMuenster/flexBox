classdef L1InfdataTerm < L1InfdataTermOperator
    methods
        function obj = L1InfdataTerm(alpha,f)
            obj = obj@L1InfdataTermOperator(alpha,identityOperator(numel(f)),f);
        end
    end
end
