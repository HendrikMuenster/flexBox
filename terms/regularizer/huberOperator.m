%represents the term
%\alpha |A u|_{H_\epsilon} (huber norm)
%corresponds to one primal variable u and an arbitrary linear operator A
classdef huberOperator < basicDualizedOperator & HuberProxDual
    properties
        epsi
    end

    methods
        function obj = huberOperator(alpha,numPrimals,A,epsi,varargin)
            obj = obj@basicDualizedOperator(alpha,numPrimals,A,varargin);
            obj.epsi = epsi;
        end

    end
end
