%
classdef tildeMultiOperatorMultiDual < handle
    methods
        function yTilde(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                main.yTilde{dualNumbers(i)} = main.y{dualNumbers(i)};
            end
            
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.yTilde{dualNumbers(i)} = main.yTilde{dualNumbers(i)} + main.params.sigma{dualNumbers(i)} * (obj.operator{operatorNumber}* main.xBar{primalNumbers(j)});
                end
            end
        end
        
        function xTilde(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.xTilde{primalNumbers(j)} = main.xTilde{primalNumbers(j)} - main.params.tau{primalNumbers(j)}*(obj.operatorT{operatorNumber} * main.y{dualNumbers(i)});
                end
            end
        end
        
        function yError(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.yError{dualNumbers(i)} = main.yError{dualNumbers(i)} - obj.operator{operatorNumber}* (main.x{primalNumbers(j)}-main.xOld{primalNumbers(j)});
                end
            end
        end
        
        function xError(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:numel(primalNumbers)
                    operatorNumber = numel(primalNumbers)*(i-1) + j;
                    main.xError{primalNumbers(j)} = main.xError{primalNumbers(j)} - obj.operatorT{operatorNumber} * (main.y{dualNumbers(i)} - main.yOld{dualNumbers(i)});
                end
            end
        end
    end
end