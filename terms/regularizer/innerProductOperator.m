%represents the term
%\alpha <b,Au>
%corresponds to one primal variable u and an arbitrary linear operator A
classdef innerProductOperator < basicDualizedOperator & innerProductProxDual
    properties
        b
    end

    methods
        function obj = innerProductOperator(alpha,A,b,varargin)
            obj = obj@basicDualizedOperator(alpha,1,A,varargin);

            if (numel(b) == size(A,1))
                obj.b{1} = b(:);

            else
                error(['Input vector b should have length ',size(A,1),' which is size(A,1)']);
            end
        end
    end
end
