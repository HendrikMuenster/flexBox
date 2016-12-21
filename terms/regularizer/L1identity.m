%represents the term
%\alpha |Iu|_1,1 where I is the identity operator
%in the L_1,1 norm
%corresponds to one primal variable u
classdef L1identity < basicIdentity & L1AnisoProxDual
    properties
    end

    methods
        function obj = L1identity(alpha,dims)
            obj = obj@basicIdentity(alpha,dims);
        end

    end
end
