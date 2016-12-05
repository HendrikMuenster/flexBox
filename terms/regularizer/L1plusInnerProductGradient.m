%represents the term
%\alpha ( \|\nabla u\|_1,1 + <b,\nabla u>)
%corresponds to one primal variable u 
classdef L1plusInnerProductGradient < basicGradient & L1IsoProxDualShift
    properties
        b
        r
    end

    methods
        function obj = L1plusInnerProductGradient(alpha,dims,b,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);

            if (numel(b) == numel(dims)*prod(dims))
                for i=1:numel(dims)
                    bTmp = b( (i-1)*prod(dims)+1 : i*prod(dims) );

                    obj.b{i} = bTmp(:);
                    obj.r{i} = bTmp(:);
                end
            else
                error(['Input vector b should have length ',num2str(numel(dims)*prod(dims)),' which is numel(dims) * prod(dims)']);
            end
        end
    end
end
