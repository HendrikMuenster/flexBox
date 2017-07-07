classdef LInfdataTerm < LInfdataTermOperator
    methods
        function obj = LInfdataTerm(alpha,f)
            obj = obj@LInfdataTermOperator(alpha,identityOperator(numel(f)),f);
        end
    end
end
