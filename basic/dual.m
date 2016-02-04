classdef dual < handle
    properties
        listDual
        correspondingPrimal
        numDual
    end
    
    methods
        function obj = dual
            obj.listDual = {};
            obj.correspondingPrimal = {};
            obj.numDual = 0;
        end
        
        function add(obj,dualPart,varargin)
            if (nargin<3)
                varargin{1} = 1; %then this dual variable corresponds to the first primal
            end
            
            obj.numDual = obj.numDual + 1;
            
            obj.listDual{obj.numDual} = dualPart;
            obj.correspondingPrimal{obj.numDual} = varargin{1};
        end
        
        function prox(obj,p,sigma)
            for j=1:numel(obj.listDual)
                obj.listDual{j}.applyProx(p,sigma,obj.correspondingPrimal{j});
            end
        end
    end
    
    
    
end
