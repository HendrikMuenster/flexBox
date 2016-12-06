%represents the term \alpha <Au,b> for one primal variable u
classdef innerProductOperator < basicDualizedOperator & innerProductProxDual
    properties
        b
    end
    
    methods
        function obj = innerProductOperator(alpha,A,b,varargin)
            if (iscell(A))
                numPrimals = numel(A);
                
                if (numel(b) ~= size(A{1},1))
                    error(['Input vector b should have length ',size(A{1},1),' which is size(A{1},1)']);
                end
            else
                numPrimals = 1;
                if (numel(b) ~= size(A,1))
                    error(['Input vector b should have length ',size(A,1),' which is size(A,1)']);
                end
            end
            
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
            
            obj.b{1} = b(:);
        end
    end
end