classdef L1AnisoProxDualShift < basicProx
    properties (Access = public, Abstract = true)
    end
    
    methods
        function obj = L1AnisoProxDual(varargin)
            
        end

        function applyProx(obj,main,dualNumbers,~)
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumbers(i)}));
            end
        end
        
    end
end