function [ result ] = segmentation(image,numberOfLabels,alpha,regularizer,varargin )
    vararginParser;
    dimsImage = size(image);
    
    % additional optional input:
    % add 'lastDim','feature' or 'lastDim','data' to decide how to 
    % treat the last dimension; if set to 'feature', then the mean 
    % values be a matrix with size(meanValues,2) = size(f,lastDim)
    % by default, lastDim is set to data
    if (~exist('lastDim','var'))
        lastDim = 'data';
    end
    
    if (strcmp(lastDim,'feature'))
        labels = rand(numberOfLabels,dimsImage(end));
        dims = dimsImage(1:end-1);
    else %lastDim is data
        labels = rand(numberOfLabels,1);
        dims = dimsImage;
    end
    
    
    %init data term
    data = labelingTerm(dims,1,image,labels,'lastDim',lastDim);
    
    %init regularizer
    if (strcmp(regularizer,'L2'))
        regularizer = L2gradient(alpha,dims);
    elseif(strcmp(regularizer,'AnisoTV'))
        regularizer = L1gradientAniso(alpha,dims);
    elseif(strcmp(regularizer,'IsoTV'))
        regularizer = L1gradientIso(alpha,dims);
    elseif(strcmp(regularizer,'HuberTV'))
        if (nargin < 6)
            error('No weight for the Huber norm given!');
        else
            epsi = varargin{1};
        end
        regularizer = L1gradientHuber(alpha,dims,epsi);
    else
        error('Please choose valid regularizer');
    end
    
    %init henPDT
    main = flexbox;
    for i=1:numberOfLabels
        main.addPrimalVar(size(image));
    end
    main.addTerm(data,1:numberOfLabels);
    
    %add regularizer for all labels
    for i=1:numberOfLabels
        main.addTerm(regularizer,i);
    end
    
    main.params.maxIt = 250; %do only 250 it
    
    err = 1;secondError = 1;
    while (err+secondError) > 1e-5
        %save old labels
        labelsOld = labels;
        
        %do some iterations
        main.runAlgorithm;
        
        %update labels
        imageShow = zeros([dimsImage(1:2),3]);
        for i=1:numberOfLabels
            if (strcmp(lastDim,'feature'))
                otherdims = repmat({':'},1,ndims(image)-1);
                for j=1:dimsImage(end)
                    if (dimsImage(end)==3) %RGB image
                        imageShow(:,:,j) = imageShow(:,:,j) + reshape(main.x{i},dims)*labels(i,j);
                    end
                    imageTmp = image(otherdims{:},j);

                    set = imageTmp(main.x{i} > 0.5);
                    if (numel(set)>0)
                        labels(i,j) = mean(set(:));
                    end
                end
            else
                if (ndims(image) == 2)%normal image
                    %create image to show
                    for j=1:3
                        imageShow(:,:,j) = imageShow(:,:,j) + reshape(main.x{i},dims)*rand(1);
                    end
                end
                
                set = image(main.x{i} > 0.5);
                if (numel(set)>0)
                    labels(i) = mean(set(:));
                end
            end
        end
        
        figure(456);imagesc(imageShow);axis image;drawnow;
        
        %update error
        err = (labels - labelsOld).^2;
        err = sum(err(:)) / numel(err);
        
        secondError = main.calculateError;
        
    end
    
    if (strcmp(lastDim,'feature'))
        otherdims = repmat({':'},1,ndims(image)-1);
        for i=1:numberOfLabels
            result(otherdims{:},i) = reshape(main.x{i},dims);
        end
    else
        otherdims = repmat({':'},1,ndims(image));
        for i=1:numberOfLabels
            result(otherdims{:},i) = reshape(main.x{i},dims);
        end
    end
    
end

