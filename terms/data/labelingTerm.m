%labeling term
classdef labelingTerm < primalPart
    
    properties
        f
        numLabels
        weights
        tmpList
        lastDim
    end

    methods
        function obj = labelingTerm(alpha,f,labels,varargin)
            vararginParser;
            %number of primal variables to be created

            if length(labels) == numel(labels)
                numL = numel(labels);
            else
                numL = size(labels,1);
            end
            
            obj = obj@primalPart(alpha);%number of primal variables equal to the number ob labels
            
            obj.numPrimals = numL;
            
            %save input image f
            obj.f = f;
            %save number of possible lables
            obj.numLabels = numL;
            
            % additional optional input:
            % add 'lastDim','feature' or 'lastDim','data' to decide how to 
            % treat the last dimension; if set to 'feature', then the mean 
            % values be a matrix with size(meanValues,2) = size(f,lastDim)
            % by default, lastDim is set to data
            if (~exist('lastDim','var'))
                obj.lastDim = 'data';
            else
                obj.lastDim = lastDim;
            end
            
            %assigns labels
            obj.setLabels(labels);
            
            %array used for simplex projection algorithm
            obj.tmpList = zeros(numel(obj.weights{1}),numL);
        end
        
        function setLabels(obj,labels)
            %decide whether last dimension stands for RGB or not
            if (strcmp(obj.lastDim,'feature'))
                otherdims = repmat({':'},1,ndims(obj.f)-1);
                
                for i=1:size(labels,1) %number of rows equals number of labels
                    obj.weights{i} = 0;
                    for j=1:size(labels,2) %number of cols equals number of features
                        obj.weights{i} = obj.weights{i} + (obj.f(otherdims{:}, j)-labels(i,j)).^2;
                    end
                    obj.weights{i} = sqrt(obj.weights{i});
                    obj.weights{i} = obj.weights{i}(:);
                end
            elseif (strcmp(obj.lastDim,'data'))
                for i=1:numel(labels)
                    obj.weights{i} = (obj.f-labels(i)).^2;
                    obj.weights{i} = obj.weights{i}(:);
                end
            end
        end
        
        function init(obj,myNumber,main)

        end
           
        function applyProx(obj,main,primalNumber)
            % (tildeX - tau*weights) should be projected onto simplex
            %
            % Uses the following algorithm:
            % Projection Onto A Simplex by Yunmei Chen and Xiaojing Ye, 2011
            % http://arxiv.org/pdf/1101.6081v2.pdf
            %
            
            for i=1:numel(primalNumber)
                obj.tmpList(:,i) = (main.xTilde{primalNumber(i)} - main.params.tau{primalNumber(i)} * obj.weights{i});
            end
            
            %do projection
            ySort = sort(obj.tmpList,2,'ascend');

            tMaxValues = zeros(size(obj.tmpList,1),1);

            for i=obj.numLabels-1:-1:1
                tMaxTmp = (sum(ySort(:,i+1:end),2)-1) / (obj.numLabels-i);

                tMaxTmp(tMaxValues~=0) = -Inf;

                trueValues = (tMaxTmp >= ySort(:,i));
                if (sum(trueValues)>0)
                    tMaxValues( trueValues) = tMaxTmp(trueValues);
                end
            end

            tMaxTmp = (sum(ySort(:,1:end),2)-1) / obj.numLabels;

            tMaxValues(tMaxValues==0) = tMaxTmp(tMaxValues==0);
            
            zVec = zeros(size(obj.tmpList,1),1);
            %update primal variables
            for i=1:numel(primalNumber)
                main.x{primalNumber(i)} = max([obj.tmpList(:,i)-tMaxValues,zVec],[],2);
            end
        end
    end
end