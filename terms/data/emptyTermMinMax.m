%prox for G = alpha  |u-f|
classdef emptyTermMinMax < primalPart & MinMaxProxPrimal

    methods
        function obj = emptyTermMinMax(minVal,maxVal)
            obj = obj@primalPart(1);%only one primal variable
            
            obj.minVal = minVal;
            obj.maxVal = maxVal;
        end
    end
end