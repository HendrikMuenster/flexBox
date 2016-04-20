function [ operator] = generateDownsamplingNew( targetDimensions,factor )


    nPxTarget = prod(targetDimensions);

    mask = reshape(1:nPxTarget,targetDimensions);


    mask = kron(mask,ones(factor));

    weight = 1/factor^2;

    operator = sparse(mask(:),1:(nPxTarget*factor*factor),weight);
end