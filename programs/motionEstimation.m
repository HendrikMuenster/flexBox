function [ result ] = motionEstimation(image1,image2,alpha,dataterm,regularizer,varargin )
    dimsU = size(image1);
    
    %init data term
    if (strcmp(dataterm,'L2'))
        data = L2opticalFlowTerm(image1,image2,1);
    elseif(strcmp(dataterm,'L1'))
        data = L1opticalFlowTerm(image1,image2,1);
    else
        error('Please choose valid data term');
    end
    
    %init regularizer
    if (strcmp(regularizer,'L2'))
        regularizer = L2gradient(alpha,dimsU);
    elseif(strcmp(regularizer,'AnisoTV'))
        regularizer = AnisoTV(alpha,dimsU);
    elseif(strcmp(regularizer,'IsoTV'))
        regularizer = IsoTV(alpha,dimsU);
    elseif(strcmp(regularizer,'HuberTV'))
        if (nargin < 6)
            error('No weight for the Huber norm given!');
        else
            epsi = varargin{1};
        end
        regularizer = HuberTV(alpha,dimsU,epsi);
    else
        error('Please choose valid regularizer');
    end
    
    %init henPDT
    main = flexbox;
    main.addPrimal(data);
    
    %add regularizer for both components
    main.addDual(regularizer,1);
    main.addDual(regularizer,2);
    
    main.runAlgorithm;
    
    v1 = reshape(main.x{1},dimsU);
    v2 = reshape(main.x{2},dimsU);
    result = cat(3,v1,v2);
end

