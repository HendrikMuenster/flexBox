function A  = warpingOperator(dimsU,flow,varargin)

if (nargin < 3)
    method = 1;
else
    method = varargin{1};
end

[meshX,meshY] = meshgrid(1:dimsU(2),1:dimsU(1));

indexList = sub2ind(dimsU, meshY, meshX);

indexList = reshape(indexList,dimsU);

targetPoint2Mat = meshX + flow(:,:,1);
targetPoint1Mat = meshY + flow(:,:,2);

x1 = floor(targetPoint1Mat) - 1;
x2 = x1 + 1;
x3 = x2 + 1;
x4 = x3 + 1;

y1 = floor(targetPoint2Mat) - 1;
y2 = y1 + 1;
y3 = y2 + 1;
y4 = y3 + 1;

v2 = targetPoint1Mat - x2;
v1 = targetPoint2Mat - y2;

indicator = ones(dimsU);
indicator(x1 < 1) = 0;
indicator(y1 < 1) = 0;
indicator(x4 > dimsU(1)) = 0;
indicator(y4 > dimsU(2)) = 0;

indicator = indicator > 0;

v1 = v1(indicator);
v2 = v2(indicator);

if (method == 1)

    listJ2 = [indexList(sub2ind(dimsU, x1(indicator(:)), y1(indicator(:))))];
    %listVal2 = (v2.^3/4 - v2.^2/2 + v2/4).*v1.^3 + (- v2.^3/2 + v2.^2 - v2/2).*v1.^2 + (v2.^3/4 - v2.^2/2 + v2/4).*v1;
    listVal2 = (v1.*v2.*(v1 - 1).^2.*(v2 - 1).^2)/4;

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x2(indicator(:)), y1(indicator(:))))];
    %listVal2 = [listVal2,(- (3.*v2.^3)/4 + (5.*v2.^2)/4 - 1/2).*v1.^3 + ((3*v2.^3)/2 - (5*v2.^2)/2 + 1).*v1.^2 + (- (3*v2.^3)/4 + (5*v2.^2)/4 - 1/2).*v1];
    listVal2 = [listVal2,-(v1.*(v1 - 1).^2.*(3.*v2.^3 - 5.*v2.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x3(indicator(:)), y1(indicator(:))))];
    %listVal2 = [listVal2,((3*v2.^3)/4 - v2.^2 - v2/4).*v1.^3 + (- (3*v2.^3)/2 + 2*v2.^2 + v2/2).*v1.^2 + ((3*v2.^3)/4 - v2.^2 - v2/4).*v1];
    listVal2 = [listVal2,-(v1.*v2.*(v1 - 1).^2.*(- 3.*v2.^2 + 4.*v2 + 1))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x4(indicator(:)), y1(indicator(:))))];
    %listVal2 = [listVal2,(v2.^2/4 - v2.^3/4).*v1.^3 + (v2.^3/2 - v2.^2/2).*v1.^2 + (v2.^2/4 - v2.^3/4).*v1];
    listVal2 = [listVal2,-(v1.*v2.^2.*(v1 - 1).^2.*(v2 - 1))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x1(indicator(:)), y2(indicator(:))))];
    %listVal2 = [listVal2,v1.^2.*((5.*v2.^3)/4 - (5.*v2.^2)/2 + (5.*v2)/4) - v1.^3.*((3.*v2.^3)/4 - (3.*v2.^2)/2 + (3.*v2)/4) - v2/2 + v2.^2 - v2.^3/2];
    listVal2 = [listVal2,-(v2.*(v2 - 1).^2.*(3.*v1.^3 - 5.*v1.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x2(indicator(:)), y2(indicator(:))))];
    %listVal2 = [listVal2,(3.*v2.^3)/2 - (5.*v2.^2)/2 + v1.^3.*((9.*v2.^3)/4 - (15.*v2.^2)/4 + 3/2) - v1.^2.*((15*v2.^3)/4 - (25.*v2.^2)/4 + 5/2) + 1];
    listVal2 = [listVal2,((3.*v1.^3 - 5.*v1.^2 + 2).*(3.*v2.^3 - 5.*v2.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x3(indicator(:)), y2(indicator(:))))];
    %listVal2 = [listVal2,v2/2 + v1.^3.*(- (9.*v2.^3)/4 + 3.*v2.^2 + (3.*v2)/4) - v1.^2.*(- (15.*v2.^3)/4 + 5.*v2.^2 + (5.*v2)/4) + 2.*v2.^2 - (3.*v2.^3)/2];
    listVal2 = [listVal2,(v2.*(- 3.*v2.^2 + 4.*v2 + 1).*(3.*v1.^3 - 5.*v1.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x4(indicator(:)), y2(indicator(:))))];
    %listVal2 = [listVal2,v2.^3/2 - v2.^2/2 - v1.^3.*(- (3.*v2.^3)/4 + (3.*v2.^2)/4) + v1.^2.*(- (5.*v2.^3)/4 + (5.*v2.^2)/4)];
    listVal2 = [listVal2,(v2.^2.*(v2 - 1).*(3.*v1.^3 - 5.*v1.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x1(indicator(:)), y3(indicator(:))))];
    %listVal2 = [listVal2,((3.*v2.^3)/4 - (3.*v2.^2)/2 + (3.*v2)/4).*v1.^3 + (- v2.^3 + 2*v2.^2 - v2).*v1.^2 + (- v2.^3/4 + v2.^2/2 - v2/4).*v1];
    listVal2 = [listVal2,-(v1.*v2.*(v2 - 1).^2.*(- 3.*v1.^2 + 4.*v1 + 1))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x2(indicator(:)), y3(indicator(:))))];
    %listVal2 = [listVal2,(- (9*v2.^3)/4 + (15*v2.^2)/4 - 3/2).*v1.^3 + (3*v2.^3 - 5*v2.^2 + 2).*v1.^2 + ((3*v2.^3)/4 - (5*v2.^2)/4 + 1/2).*v1];
    listVal2 = [listVal2,(v1.*(- 3.*v1.^2 + 4.*v1 + 1).*(3.*v2.^3 - 5.*v2.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x3(indicator(:)), y3(indicator(:))))];
    %listVal2 = [listVal2,v1.*(- (3*v2.^3)/4 + v2.^2 + v2/4) + v1.^2.*(- 3*v2.^3 + 4*v2.^2 + v2) - v1.^3.*(- (9*v2.^3)/4 + 3*v2.^2 + (3*v2)/4)];
    listVal2 = [listVal2,(v1.*v2.*(- 3.*v1.^2 + 4.*v1 + 1).*(- 3.*v2.^2 + 4.*v2 + 1))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x4(indicator(:)), y3(indicator(:))))];
    %listVal2 = [listVal2,((3*v2.^2)/4 - (3*v2.^3)/4).*v1.^3 + (v2.^3 - v2.^2).*v1.^2 + (v2.^3/4 - v2.^2/4).*v1];
    listVal2 = [listVal2,(v1.*v2.^2.*(v2 - 1).*(- 3.*v1.^2 + 4.*v1 + 1))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x1(indicator(:)), y4(indicator(:))))];
    %listVal2 = [listVal2,(- v2.^3/4 + v2.^2/2 - v2/4).*v1.^3 + (v2.^3/4 - v2.^2/2 + v2/4).*v1.^2];
    listVal2 = [listVal2,-(v1.^2.*v2.*(v1 - 1).*(v2 - 1).^2)./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x2(indicator(:)), y4(indicator(:))))];
    %listVal2 = [listVal2,((3*v2.^3)/4 - (5*v2.^2)/4 + 1/2).*v1.^3 + (- (3*v2.^3)/4 + (5*v2.^2)/4 - 1/2).*v1.^2];
    listVal2 = [listVal2,(v1.^2.*(v1 - 1).*(3.*v2.^3 - 5.*v2.^2 + 2))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x3(indicator(:)), y4(indicator(:))))];
    %listVal2 = [listVal2,v1.^3.*(- (3*v2.^3)/4 + v2.^2 + v2/4) - v1.^2.*(- (3*v2.^3)/4 + v2.^2 + v2/4)];
    listVal2 = [listVal2,(v1.^2.*v2.*(v1 - 1).*(- 3.*v2.^2 + 4.*v2 + 1))./4];

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x4(indicator(:)), y4(indicator(:))))];
    %listVal2 = [listVal2,(v2.^3/4 - v2.^2/4).*v1.^3 + (v2.^2/4 - v2.^3/4).*v1.^2];
    listVal2 = [listVal2,(v1.^2.*v2.^2.*(v1 - 1).*(v2 - 1))./4];

    listI2 = indexList(indicator);
    listI2 = repmat(listI2,1,16);
    
elseif (method == 2)

    listJ2 = [indexList(sub2ind(dimsU, x2(indicator(:)), y2(indicator(:))))];
    listVal2 = (1-v1) .* (1-v2);

    listJ2 = [listJ2,indexList(sub2ind(dimsU, x3(indicator(:)), y2(indicator(:))))];
    listVal2 = [listVal2,v2 .* (1-v1)];
    
    listJ2 = [listJ2,indexList(sub2ind(dimsU, x2(indicator(:)), y3(indicator(:))))];
    listVal2 = [listVal2,v1 .* (1-v2)];
    
    listJ2 = [listJ2,indexList(sub2ind(dimsU, x3(indicator(:)), y3(indicator(:))))];
    listVal2 = [listVal2,v1.*v2];

    listI2 = indexList(indicator);
    listI2 = repmat(listI2,1,4);
end


A = sparse(listI2(:),listJ2(:),listVal2(:));
A(1,prod(dimsU)) = 0;
A(prod(dimsU),1) = 0;

end