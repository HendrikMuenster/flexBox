%represents the term
%\alpha |Iu|_2^2 where I is the identity operator
%in the L_2^2 norm
%corresponds to one primal variable u
classdef L2identity < basicIdentity & L2proxDual
    properties
    end

    methods
        function obj = L2identity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
        end

    end
end
