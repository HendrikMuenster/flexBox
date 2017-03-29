classdef labelingProxDual < basicProx
    properties
    end
    
    methods
        function obj = labelingProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            % % (tildeX - tau*weights) should be projected onto simplex
            % %
            % % Uses the following algorithm:
            % % Projection Onto A Simplex by Yunmei Chen and Xiaojing Ye, 2011
            % % http://arxiv.org/pdf/1101.6081v2.pdf
            % %
            
			%use Moreau decomposition
			%y = tildeY - sigma * prox_{1/sigma}F(tildeY)
		
			%(1) calculate prox_{1/sigma}F(tildeY/sigma)
			
            for i=1:obj.numVars
                tmpList(:,i) = (main.yTilde{dualNumbers(i)} - obj.f{i}) ./ main.params.sigma{dualNumbers(i)};
            end
            
            %do projection
            ySort = sort(tmpList,2,'ascend');

            tMaxValues = zeros(size(tmpList,1),1);

            for i=obj.numVars-1:-1:1
                tMaxTmp = (sum(ySort(:,i+1:end),2)-1) / (obj.numVars-i);

                tMaxTmp(tMaxValues~=0) = -Inf;

                trueValues = (tMaxTmp >= ySort(:,i));
                if (sum(trueValues)>0)
                    tMaxValues( trueValues) = tMaxTmp(trueValues);
                end
            end

            tMaxTmp = (sum(ySort(:,1:end),2)-1) / obj.numVars;

            tMaxValues(tMaxValues==0) = tMaxTmp(tMaxValues==0);
            
            zVec = zeros(size(tmpList,1),1);
            %update primal variables
            for i=1:obj.numVars
                main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)} - main.params.sigma{dualNumbers(i)} .* max([tmpList(:,i)-tMaxValues,zVec],[],2);
            end
			
        end
        
    end
end