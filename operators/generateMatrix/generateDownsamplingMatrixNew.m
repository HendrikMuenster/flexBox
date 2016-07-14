function [ operator] = generateDownsamplingMatrixNew( newDims,oldDims )
    

    Xgrid = (1:oldDims(2));
    XgridNew = ( (1 + (oldDims(2)/newDims(2) - 1)/2) : oldDims(2)/newDims(2) : oldDims(2) ) ;

    Xgrid = Xgrid / 4;
    XgridNew = XgridNew / 4;
    
    
    mat1 = sparse(newDims(2),oldDims(2));
    for i=1:oldDims(2)
        vector = zeros([1,oldDims(2)]);
        vector(i) = 1;
        mat1(:,i) = interp1(Xgrid,vector,XgridNew,'pchip');
    end

    Ygrid = 1:oldDims(1);
    YgridNew = (1 + (oldDims(1)/newDims(1) - 1)/2) : oldDims(1)/newDims(1) : oldDims(1);
    
    Ygrid = Ygrid / 4;
    YgridNew = YgridNew / 4;
    
    mat2 = sparse(newDims(1),oldDims(1));
    for i=1:oldDims(1)
        vector = zeros([1,oldDims(1)]);
        vector(i) = 1;
        mat2(:,i) = interp1(Ygrid,vector,YgridNew,'pchip');
    end

     operator = kron(mat1,mat2);
     
%      [X,Y] = meshgrid(Xgrid,Ygrid);
%      [Xq,Yq] = meshgrid(XgridNew,YgridNew);
%      
%      operator2 = sparse(size(operator,1),size(operator,2));
%     for i=1:oldDims(1) * oldDims(2)
% 
%         V = zeros(size(X));
%         V(i) = 1;
%         Vq = interp2(X,Y,V,Xq,Yq,'cubic');
% 
%         operator2(:,i) = Vq(:);
% 
%     end
% 
%     operator = operator2;
     
     
     
end