%
classdef tildeMultiOperatorSingleDual < handle
    methods
        function yTilde(obj,main,dualNumbers,primalNumber)
            for i=1:obj.numVars
                main.yTilde{dualNumbers(i)} = main.y{dualNumbers(i)} + main.params.sigma{dualNumbers(i)} * (obj.operator{i}* main.xBar{primalNumber});
            end
        end
        
        function xTilde(obj,main,dualNumbers,primalNumber)
            for i=1:obj.numVars
                main.xTilde{primalNumber} = main.xTilde{primalNumber} - main.params.tau {primalNumber} * (obj.operator{i}' * main.y{dualNumbers(i)});
            end
        end
        
        function yError(obj,main,dualNumbers,primalNumber)
            for i=1:obj.numVars
                main.yError{dualNumbers(i)} = main.yError{dualNumbers(i)} - obj.operator{i} * (main.x{primalNumber}-main.xOld{primalNumber});
            end
        end
        
        function xError(obj,main,dualNumbers,primalNumber)
            for i=1:obj.numVars
                main.xError{primalNumber} = main.xError{primalNumber} - obj.operator{i}' * (main.y{dualNumbers(i)}-main.yOld{dualNumbers(i)});
            end
        end
    end
end