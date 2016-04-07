function [ operator] = generateDownsamplingMatrix( newDims,oldDims )
    
    [Xgrid] = 1:(newDims(2)-1)/(oldDims(2)-1):newDims(2);
    [XgridNew] = 1:newDims(2);

    mat1 = sparse(newDims(2),oldDims(2));
    for i=1:oldDims(2)
        vector = zeros([1,oldDims(2)]);
        vector(i) = 1;
        mat1(:,i) = interp1(Xgrid,vector,XgridNew,'PCHIP');
    end

    [Ygrid] =1:(newDims(1)-1)/(oldDims(1)-1):newDims(1);
    [YgridNew] = 1:newDims(1);

    mat2 = sparse(newDims(1),oldDims(1));
    for i=1:oldDims(1)
        vector = zeros([1,oldDims(1)]);
        vector(i) = 1;
        mat2(:,i) = interp1(Ygrid,vector,YgridNew,'PCHIP');
    end

     operator = kron(mat2,mat1);
end