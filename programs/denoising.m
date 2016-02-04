function [ result ] = denoising( f,alpha,dataterm,regularizer,varargin )
    vararginParser;
    
    dimsU = size(f);
    
    main = flexbox;
    
    %init data term
    if (strcmp(dataterm,'L2'))
        main.addPrimal(L2dataTerm(dimsU,1,f));
    elseif(strcmp(dataterm,'L1'))
        main.addPrimal(L1dataTerm(dimsU,1,f));
    elseif(strcmp(dataterm,'KL'))
        main.addPrimal(emptyTermMinMax(dimsU,0,1));
        main.addDual(KullbackLeiblerDualizedDataTerm(1,speye(prod(dimsU)),f),1);
    else
        error('Please choose valid data term');
    end
    
    
    varsDual = 1;%corresponding primal var for dual
    %init regularizer
    if (strcmp(regularizer,'L2'))
        regularizer1 = L2gradient(alpha,dimsU);
    elseif(strcmp(regularizer,'AnisoTV'))
        regularizer1 = L1gradientAniso(alpha,dimsU);
    elseif(strcmp(regularizer,'IsoTV'))
        regularizer1 = L1gradientIso(alpha,dimsU);
    elseif(strcmp(regularizer,'HuberTV'))
        if (~exist('huberWeight','var'))
            error('No weight for the Huber norm given!');
        else
            epsi = huberWeight;
        end
        regularizer1 = HuberTV(alpha,dimsU,epsi);
    elseif(strcmp(regularizer,'TV2'))
        if (~exist('alpha2','var'))
            error('No weight for second regularizer given!');
        end
            
        %add additional primal variables
        main.addPrimal(emptyDataTerm(dimsU));
        main.addPrimal(emptyDataTerm(dimsU));

        regularizer1 = IsoTVSecondOrder(alpha2,dimsU);
        regularizer2 = AnisoTV(1,dimsU);
        
        main.addDual(regularizer2,2);
        main.addDual(regularizer2,3);
        varsDual = [1,2,3];
    else
        error('Please choose valid regularizer');
    end
    
    main.addDual(regularizer1,varsDual);
    
    main.params.showPrimals = 100;
    
    %init henPDT
    main.runAlgorithm;
    
    
    
    result = main.getPrimal(1);
end

