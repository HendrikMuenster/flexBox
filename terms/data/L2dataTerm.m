% L2dataTerm(dims,alpha,f) represents an L2 data term 
% alpha/2 ||u-f||_2^2
classdef L2dataTerm < L2dataTermOperator
    methods
        function obj = L2dataTerm(alpha,f)
            obj = obj@L2dataTermOperator(alpha,identityOperator(numel(f)),f);
        end
    end
end