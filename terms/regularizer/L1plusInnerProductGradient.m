%represents the term \alpha ( \|\nabla u\|_1 + <b,\nabla u>) for one primal variable u
classdef L1plusInnerProductGradient < basicGradient & L1IsoProxDualShift
    properties
        b
    end
    
    methods
        function obj = L1plusInnerProductGradient(alpha,dims,b,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
            
            if (numel(b) == numel(dims)*prod(dims))
                for i=1:numel(dims)
                    obj.b{i} = b( (i-1)*prod(dims)+1 : i*prod(dims) );
                end
            else
                error(['Input vector b should have length ',num2str(numel(dims)*prod(dims)),' which is numel(dims) * prod(dims)']);
            end
        end
    end
end