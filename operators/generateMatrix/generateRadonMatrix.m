% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.0
% Date: 2017-01-12

% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and University of Muenster not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
% 
%
%
% This function calculates the sparse matrix representation of MATLABs radon
% function. This is needed, since iradon is not exactly the discrete
% adjoint to radon, but a discretization of the analytic adjoint. 
%
% input: size of the input image and list of angles
% output: sparse matrix representation of the radon function and size array of
% the corresponding sinogram
function [ radonMatrix, sizeSinogram] = generateRadonMatrix ( imageSize,angles )
    largeNumber = 1e4;

    N = prod(imageSize);

    emptyImage = zeros(imageSize);

    Gi = zeros(largeNumber,1);
    Gj = zeros(largeNumber,1);
    Gx = zeros(largeNumber,1);

    start = 1;

    reverseStr = '';
    for i=1:N
        emptyImage(i) = 1;

        radonResult = radon(emptyImage, angles);

        listNonzeroEntries = find(radonResult(:));

        if (start + numel(listNonzeroEntries) > numel(Gi))
            Gi = [Gi;zeros(largeNumber,1)];
            Gj = [Gj;zeros(largeNumber,1)];
            Gx = [Gx;zeros(largeNumber,1)];
        end

        Gj(start:start + numel(listNonzeroEntries) - 1) = i*ones(numel(listNonzeroEntries),1);
        Gi(start:start + numel(listNonzeroEntries) - 1) = listNonzeroEntries;
        Gx(start:start + numel(listNonzeroEntries) - 1) = radonResult(listNonzeroEntries);

        emptyImage(i) = 0;

        start = start + numel(listNonzeroEntries);

        if (mod(i,1000) == 1)
            [ reverseStr ] = printToCmd( reverseStr,sprintf('Generating matrix equivalent of MATLABs Radon transform: %0.0f percent',round(100*i/N)));
        end
    end

    listZeroEntries = find(Gi==0);

    Gi(listZeroEntries) = [];
    Gj(listZeroEntries) = [];
    Gx(listZeroEntries) = [];


    radonMatrix = sparse(Gi,Gj,Gx,numel(radonResult),numel(emptyImage));
    sizeSinogram = size(radonResult);

    printToCmd( reverseStr,'');

end