%labeling term
classdef labelingTerm < basicDualizedDataterm & labelingProxDual
    
    properties
    end

    methods
        function obj = labelingTerm(alpha,f,labels,varargin)
            vararginParser;

            % additional optional input:
            % add 'lastDim','feature' or 'lastDim','data' to decide how to 
            % treat the last dimension; if set to 'feature', then the mean 
            % values be a matrix with size(meanValues,2) = size(f,lastDim)
            % by default, lastDim is set to data
			if (~exist('lastDim','var'))
                lastDim = 'data';
			end

            if length(labels) == numel(labels)
                numL = numel(labels);
            else
                numL = size(labels,1);
            end
            
            if (strcmp(lastDim,'feature'))
                sizeOperator = numel(f(:,:,1));
            else
                sizeOperator = numel(f);
            end
			
            %operator is a diagonal block matrix with identity operators
			for i=1:numL
                for j=1:numL
                    operatorNumber = (i-1)*numL + j;
                    
                    if (i==j)
                        A{operatorNumber} = identityOperator(sizeOperator);
                    else
                        A{operatorNumber} = zeroOperator(sizeOperator);
                    end
                end
			end
			
			weight = labelingTerm.calculateWeight(labels,f,lastDim);
			
			obj = obj@basicDualizedDataterm(alpha,numL,A,weight);

            %array used for simplex projection algorithm
            %obj.tmpList = zeros(numel(obj.weights{1}),numL);
        end
    end
    methods(Static)
        function weights = calculateWeight(labels,f,lastDim)
            %decide whether last dimension stands for RGB or not
            if (strcmp(lastDim,'feature'))
                otherdims = repmat({':'},1,ndims(f)-1);
                
                for i=1:size(labels,1) %number of rows equals number of labels
                    weights{i} = 0;
                    for j=1:size(labels,2) %number of cols equals number of features
                        weights{i} = weights{i} + (f(otherdims{:}, j)-labels(i,j)).^2;
                    end
                    weights{i} = sqrt(weights{i});
                    weights{i} = weights{i}(:);
                end
            elseif (strcmp(lastDim,'data'))
                for i=1:numel(labels)
                    weights{i} = (f-labels(i)).^2;
                    weights{i} = weights{i}(:);
                end
            end
        end
        
           
        % function applyProx(obj,main,primalNumber)
            % % (tildeX - tau*weights) should be projected onto simplex
            % %
            % % Uses the following algorithm:
            % % Projection Onto A Simplex by Yunmei Chen and Xiaojing Ye, 2011
            % % http://arxiv.org/pdf/1101.6081v2.pdf
            % %
            
            % for i=1:numel(primalNumber)
                % obj.tmpList(:,i) = (main.xTilde{primalNumber(i)} - main.params.tau{primalNumber(i)} * obj.weights{i});
            % end
            
            % %do projection
            % ySort = sort(obj.tmpList,2,'ascend');

            % tMaxValues = zeros(size(obj.tmpList,1),1);

            % for i=obj.numLabels-1:-1:1
                % tMaxTmp = (sum(ySort(:,i+1:end),2)-1) / (obj.numLabels-i);

                % tMaxTmp(tMaxValues~=0) = -Inf;

                % trueValues = (tMaxTmp >= ySort(:,i));
                % if (sum(trueValues)>0)
                    % tMaxValues( trueValues) = tMaxTmp(trueValues);
                % end
            % end

            % tMaxTmp = (sum(ySort(:,1:end),2)-1) / obj.numLabels;

            % tMaxValues(tMaxValues==0) = tMaxTmp(tMaxValues==0);
            
            % zVec = zeros(size(obj.tmpList,1),1);
            % %update primal variables
            % for i=1:numel(primalNumber)
                % main.x{primalNumber(i)} = max([obj.tmpList(:,i)-tMaxValues,zVec],[],2);
            % end
        % end
    end
end