%
classdef tildeSingleOperator < handle
    methods
        function yTilde(obj,main,dualNumber,primalNumber)
            main.yTilde{dualNumber} = main.y{dualNumber} + main.params.sigma{dualNumber} * (obj.operator{1}* main.xBar{primalNumber});
        end
        
        function xTilde(obj,main,dualNumber,primalNumber)
            main.xTilde{primalNumber} = main.xTilde{primalNumber} - main.params.tau {primalNumber} * (obj.operator{1}' * main.y{dualNumber});
        end
        
        function yError(obj,main,dualNumber,primalNumber)
            main.yError{dualNumber} = main.yError{dualNumber} - obj.operator{1} * (main.x{primalNumber}-main.xOld{primalNumber});
        end
        
        function xError(obj,main,dualNumber,primalNumber)
            main.xError{primalNumber} = main.xError{primalNumber} - obj.operator{1}' * (main.y{dualNumber}-main.yOld{dualNumber});
        end
    end
end