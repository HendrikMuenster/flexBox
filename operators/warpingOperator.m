function A  = warpingOperator(dimsU,flow,varargin)

listI = [];
listJ = [];
listVal = [];

%%
[meshX,meshY] = meshgrid(1:dimsU(2),1:dimsU(1));

indexList = sub2ind(dimsU, meshY, meshX);

indexList = reshape(indexList,dimsU);

if (nargin > 2)
targetValue = interp2(varargin{2},meshX+flow(:,:,1),meshY + flow(:,:,2));
end
%%
for i=1:dimsU(1)
    for j=1:dimsU(2)
        skip = 0;

        targetPoint1 = i + flow(i,j,2);
        targetPoint2 = j + flow(i,j,1);

        x(1) = floor(targetPoint1) - 1;
        x(2) = x(1) + 1;
        x(3) = x(2) + 1;
        x(4) = x(3) + 1;

        y(1) = floor(targetPoint2) - 1;
        y(2) = y(1) + 1;
        y(3) = y(2) + 1;
        y(4) = y(3) + 1;

        v2 = targetPoint1 - x(2);
        v1 = targetPoint2 - y(2);
        

        numRows = indexList(i,j);

        if (x(2) < 1 || y(2) < 1 || x(3) > dimsU(1) || y(3) > dimsU(2))
            skip = 1;
        elseif (x(1) < 1 || y(1) < 1 || x(4) > dimsU(1) || y(4) > dimsU(2))
            %use linear interpolation if possible

            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(2),y(2));
            listVal(end+1) = (1-v1) * (1-v2);

            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(3),y(2));
            listVal(end+1) = v2 * (1-v1);

            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(2),y(3));
            listVal(end+1) = v1 * (1-v2);

            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(3),y(3));
            listVal(end+1) = v1*v2;

            %skip the rest
            skip = 1;
        end

        if (~skip)

            weight = ones(4);
            if (nargin > 2)
                templateResult = varargin{1};
                interpolatedFunction = varargin{2};
%i
%j
                %adjust interpolation to intensity values
                for nx = 1:4
                    for ny = 1:4
                        
                        %intensityDifference(ny,nx) = exp(-abs(templateResult(i,j) - interpolatedFunction(x(nx),y(ny)))^2 / 2);
                        intensityDifference(ny,nx) = exp(-abs(templateResult(i,j) - targetValue(i,j))^2);

                    end
                end
intensityDifference = 16*intensityDifference / sum(intensityDifference(:));

%                 if (i>20 && j>20)
%                  figure(1010);imagesc(intensityDifference);colorbar;pause
%                  end
weight = intensityDifference;
            end
            
            %x = xTmp;
            
            %p00
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(1),y(1));
            listVal(end+1) = weight(1,1) * ( (v2^3/4 - v2^2/2 + v2/4)*v1^3 + (- v2^3/2 + v2^2 - v2/2)*v1^2 + (v2^3/4 - v2^2/2 + v2/4)*v1);

            %p01
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(2),y(1));
            listVal(end+1) = weight(2,1) * ( (- (3*v2^3)/4 + (5*v2^2)/4 - 1/2)*v1^3 + ((3*v2^3)/2 - (5*v2^2)/2 + 1)*v1^2 + (- (3*v2^3)/4 + (5*v2^2)/4 - 1/2)*v1);

            %p02
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(3),y(1));
            listVal(end+1) = weight(3,1) * ( ((3*v2^3)/4 - v2^2 - v2/4)*v1^3 + (- (3*v2^3)/2 + 2*v2^2 + v2/2)*v1^2 + ((3*v2^3)/4 - v2^2 - v2/4)*v1);

            %p03
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(4),y(1));
            listVal(end+1) = weight(4,1) * ( (v2^2/4 - v2^3/4)*v1^3 + (v2^3/2 - v2^2/2)*v1^2 + (v2^2/4 - v2^3/4)*v1);

            %p10
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(1),y(2));
            listVal(end+1) = weight(1,2) * ( v1^2*((5*v2^3)/4 - (5*v2^2)/2 + (5*v2)/4) - v1^3*((3*v2^3)/4 - (3*v2^2)/2 + (3*v2)/4) - v2/2 + v2^2 - v2^3/2);

            %p11
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(2),y(2));
            listVal(end+1) = weight(2,2) * ( (3*v2^3)/2 - (5*v2^2)/2 + v1^3*((9*v2^3)/4 - (15*v2^2)/4 + 3/2) - v1^2*((15*v2^3)/4 - (25*v2^2)/4 + 5/2) + 1);

            %p12
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(3),y(2));
            listVal(end+1) = weight(3,2) * ( v2/2 + v1^3*(- (9*v2^3)/4 + 3*v2^2 + (3*v2)/4) - v1^2*(- (15*v2^3)/4 + 5*v2^2 + (5*v2)/4) + 2*v2^2 - (3*v2^3)/2);

            %p13
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(4),y(2));
            listVal(end+1) = weight(4,2) * ( v2^3/2 - v2^2/2 - v1^3*(- (3*v2^3)/4 + (3*v2^2)/4) + v1^2*(- (5*v2^3)/4 + (5*v2^2)/4));

            %p20
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(1),y(3));
            listVal(end+1) = weight(1,3) * ( ((3*v2^3)/4 - (3*v2^2)/2 + (3*v2)/4)*v1^3 + (- v2^3 + 2*v2^2 - v2)*v1^2 + (- v2^3/4 + v2^2/2 - v2/4)*v1);

            %p21
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(2),y(3));
            listVal(end+1) = weight(2,3) * ( (- (9*v2^3)/4 + (15*v2^2)/4 - 3/2)*v1^3 + (3*v2^3 - 5*v2^2 + 2)*v1^2 + ((3*v2^3)/4 - (5*v2^2)/4 + 1/2)*v1);

            %p22
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(3),y(3));
            listVal(end+1) = weight(3,3) * ( v1*(- (3*v2^3)/4 + v2^2 + v2/4) + v1^2*(- 3*v2^3 + 4*v2^2 + v2) - v1^3*(- (9*v2^3)/4 + 3*v2^2 + (3*v2)/4));

            %p23
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(4),y(3));
            listVal(end+1) = weight(4,3) * ( ((3*v2^2)/4 - (3*v2^3)/4)*v1^3 + (v2^3 - v2^2)*v1^2 + (v2^3/4 - v2^2/4)*v1);

            %p30
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(1),y(4));
            listVal(end+1) = weight(1,4) * ( (- v2^3/4 + v2^2/2 - v2/4)*v1^3 + (v2^3/4 - v2^2/2 + v2/4)*v1^2);

            %p31
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(2),y(4));
            listVal(end+1) = weight(2,4) * ( ((3*v2^3)/4 - (5*v2^2)/4 + 1/2)*v1^3 + (- (3*v2^3)/4 + (5*v2^2)/4 - 1/2)*v1^2);

            %p32
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(3),y(4));
            listVal(end+1) = weight(3,4) * ( v1^3*(- (3*v2^3)/4 + v2^2 + v2/4) - v1^2*(- (3*v2^3)/4 + v2^2 + v2/4));

            %p33
            listI(end+1) = numRows;
            listJ(end+1) = indexList(x(4),y(4));
            listVal(end+1) = weight(4,4) * ( (v2^3/4 - v2^2/4)*v1^3 + (v2^2/4 - v2^3/4)*v1^2) ;
        end
    end
end
%%
A = sparse(listI,listJ,listVal);
A(1,prod(dimsU)) = 0;
A(prod(dimsU),1) = 0;

end