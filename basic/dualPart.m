classdef dualPart < functionalPart
    properties
        operator  %operator
        operatorT %transposed operator
        numVars   %num of dual variables
        length    %vector of length numVars containing the length of corresponding dual variable
        myTau     %primal dual normalization factor for primal part
        mySigma   %primal dual normalization factor for dual part
    end

    methods
        function obj = dualPart(alpha)
            obj = obj@functionalPart(alpha);
        end
        
        function init(obj)
            for i=1:obj.numPrimals
                obj.myTau{i} = 0;
            end
            for i=1:numel(obj.operator) / obj.numPrimals
                obj.mySigma{i} = 0;
            end
            
            for i=1:numel(obj.operator) / obj.numPrimals %for all dual variables
                for j=1:obj.numPrimals %for all primals
                    opNum = (i-1)*obj.numPrimals + j;
                    
                    if (isa(obj.operator{opNum},'basicOperator'))
                        %this method must be implemented by every
                        %custom operator
                        obj.mySigma{i} = obj.mySigma{i} + obj.operator{opNum}.getMaxRowSumAbs();
                        obj.myTau{j} = obj.myTau{j} + obj.operatorT{opNum}.getMaxRowSumAbs();
                    else
                        obj.mySigma{i} = obj.mySigma{i} + max(sum(abs(opTmp),1));
                        obj.myTau{j} = obj.myTau{j} + max(sum(abs(opTmp),2));
                    end
                end
            end
        end
    end

    methods (Abstract)
        applyProx(obj)
    end
end
