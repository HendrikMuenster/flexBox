%represents the term
%\alpha |\nabla u|_F (frobenius norm)
%corresponds to one primal variable u
classdef frobeniusGradient < basicGradient & FrobeniusProxDual
    properties
    end

    methods
        function obj = frobeniusGradient(alpha,dims,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);

            obj.CPPsupport = 1;
        end
    end
end
