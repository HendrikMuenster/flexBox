%represents the term
%\alpha |\nabla u|_{H_\epsilon} (huber norm)
%corresponds to one primal variable u
classdef huberGradient < basicGradient & HuberProxDual
    properties
        epsi
    end

    methods
        function obj = huberGradient(alpha,dims,epsi,varargin)
            obj = obj@basicGradient(alpha,dims,varargin);
            obj.epsi = epsi;
        end

    end
end
