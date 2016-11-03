classdef functionalPart < handle

    properties
        factor      %term weight
		    CPPsupport  %flag if deriving term supports C++
        numPrimals  %corresponding number of primals
    end

    methods
        function obj = functionalPart(alpha)
            obj.factor = alpha;

            %default values
			      obj.CPPsupport = 0;
            obj.numPrimals = 1;
        end

        %function call to preallocate usefull things
        function init(varargin)
        end
    end



    methods (Abstract)
        applyProx(obj)  %abstract method for prox function (see primal dual algorithm)
    end
end
