%prox for G = alpha  |u-f|
classdef emptyDataTerm < primalPart & IdentityProxPrimal

    methods
        function obj = emptyDataTerm()
            obj = obj@primalPart(1);%only one primal variable
            
            obj.CPPsupport = 1;
        end
    end
end