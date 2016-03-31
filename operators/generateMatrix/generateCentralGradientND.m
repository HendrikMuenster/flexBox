% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.0
% Date: 2015-06-17

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
% The gradient in dimension i can be built by applying kronecker products on the sequence
% eye(dims(n)),...,eye(dims(i+1)),Dx(dims(i)),eye(dims(i-1)),...,eye(dims(1))
function [ gradientComplete] = generateCentralGradientND( dims,stepsize,varargin )
    gradientComplete = [];
    
    if (nargin > 2)
        %if third agument exists, then this specifies a specific dimension
        i = varargin{1};
        
        gradMat = 1/(2*stepsize(i)) .* spdiags([-ones(dims(i), 1),zeros(dims(i), 1),ones(dims(i), 1)], -1:1, dims(i), dims(i));
        gradMat(1,:) = 0;
        gradMat(end,:) = 0;

        if (i==1)
            grad = gradMat;
        else
            grad = speye(dims(1));
        end

        for j=2:numel(stepsize)
            if (j==i)
                grad = kron(gradMat,grad);
            else
                grad = kron(speye(dims(j)),grad);
            end
        end
        gradientComplete = [gradientComplete;grad];
    else
        for i=1:numel(stepsize)
            gradMat = 1/(2*stepsize(i)) .* spdiags([-ones(dims(i), 1),zeros(dims(i), 1),ones(dims(i), 1)], -1:1, dims(i), dims(i));
            gradMat(1,:) = 0;
            gradMat(end,:) = 0;

            if (i==1)
                grad = gradMat;
            else
                grad = speye(dims(1));
            end

            for j=2:numel(stepsize)
                if (j==i)
                    grad = kron(gradMat,grad);
                else
                    grad = kron(speye(dims(j)),grad);
                end
            end
            gradientComplete = [gradientComplete;grad];
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    

end