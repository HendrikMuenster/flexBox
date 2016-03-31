%term respresenting an indicator function on the interval [minVal,maxVal]
classdef emptyTermMinMax < primalPart & MinMaxProxPrimal

    methods
        function obj = emptyTermMinMax(minVal,maxVal)
            obj = obj@primalPart(1);%only one primal variable
            
            obj.minVal = minVal;
            obj.maxVal = maxVal;
        end
    end
end