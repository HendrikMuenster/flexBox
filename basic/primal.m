classdef primal < handle
    properties
        listPrimal
        nP %number of primal variables
        indexP1
        indexP2
    end
    
    methods
        function obj = primal
            obj.listPrimal = {};
            obj.nP = 0;
            obj.indexP1 = {};
            obj.indexP2 = {};
        end
        
        function add(obj,primalPart)
            obj.listPrimal{end + 1} = primalPart;
            
            for i=1:primalPart.number
                obj.indexP1{obj.nP + i} = numel(obj.listPrimal); %index of the primal term
                obj.indexP2{obj.nP + i} = i; %index of the inner variable
            end
            
            obj.nP = obj.nP + primalPart.number;
        end
        
        function prox(obj,d,tau)
            for j=1:numel(obj.listPrimal)
                obj.listPrimal{j}.applyProx(d,tau);
            end
        end
    end
    
    
    
end

