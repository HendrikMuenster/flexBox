% basis class for optical flow data terms u_t + \nabla u\cdot v
classdef basicOpticalFlow < basicDualizedDataterm
    properties
        thisImage1;
        thisImage2;
        termType %brightnessConstancy (default) or gradientConstancy
        constancyDimension
    end
    
    methods
        function obj = basicOpticalFlow(alpha,image1,image2,varargin)
            if (nargin > 3 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            initVar('discretization','interpolated');
            initVar('v1Tilde',zeros(size(image1)));
            initVar('v2Tilde',zeros(size(image1)));
            initVar('termType','brightnessConstancy');
            
            if (strcmp(termType,'gradientConstancy'))
                if (~exist('constancyDimension','var'))
                    error('Please specify constancy dimension by adding ''constancyDimension'',''dim''');
                end
            else
                constancyDimension = 0;
            end
            
            
            [ux,uy,ut] = basicOpticalFlow.generateDatatermParts(discretization,image1,image2,v1Tilde,v2Tilde,termType,constancyDimension);
            
            A{1} = diagonalOperator(ux);
            A{2} = diagonalOperator(uy);
            
            obj = obj@basicDualizedDataterm(alpha,2,A,-ut(:),varargin);
            
            obj.thisImage1 = image1;
            obj.thisImage2 = image2;
            obj.termType = termType;
            obj.constancyDimension = constancyDimension;
        end
        
        function warpDataterm(obj,v1Tilde,v2Tilde,varargin)
            initVar('discretization','interpolated');
            
            v1Tilde = reshape(v1Tilde,size(obj.thisImage1));
            v2Tilde = reshape(v2Tilde,size(obj.thisImage1));
            
            [ux,uy,ut] = obj.generateDatatermParts(discretization,obj.thisImage1,obj.thisImage2,v1Tilde,v2Tilde,obj.termType,obj.constancyDimension);
            
            obj.operator{1} = diagonalOperator(ux);
            obj.operatorT{1} = diagonalOperator(ux);
            
            obj.operator{2} = diagonalOperator(uy);
            obj.operatorT{2} = diagonalOperator(uy);
            
            obj.f{1} = -ut;
            
        end
        
    end
    methods(Static)
        function [ux,uy,ut] = generateDatatermParts(discretization,image1,image2,v1Tilde,v2Tilde,termType,constancyDimension)
            [I1x,I1y,I2w,I2wx,I2wy,I2wxx,I2wxy,I2wyx,I2wyy,markerOutOfGrid] = basicOpticalFlow.generateGradients(discretization,image1,image2,v1Tilde,v2Tilde);
            
            if (strcmp(termType,'brightnessConstancy'))
                ux = I2wx;
                uy = I2wy;
                ut = I2w - image1 - v1Tilde .* I2wx -  v2Tilde .* I2wy;
            elseif (strcmp(termType,'gradientConstancy'))
                if (constancyDimension == 1)
                    ux = I2wxx;
                    uy = I2wxy;
                    ut = I2wx - I1x - v1Tilde .* I2wxx -  v2Tilde .* I2wxy;
                elseif (constancyDimension == 2)
                    ux = I2wyx;
                    uy = I2wyy;
                    ut = I2wy - I1y - v1Tilde .* I2wyx -  v2Tilde .* I2wyy;
                end
            end

            if (sum(markerOutOfGrid(:)) > 0)
                ux(markerOutOfGrid) = 0;
                uy(markerOutOfGrid) = 0;
                ut(markerOutOfGrid) = 0;
            end
        end
        
        function [I1x,I1y,I2w,I2wx,I2wy,I2wxx,I2wxy,I2wyx,I2wyy,markerOutOfGrid] = generateGradients(discretization,image1,image2,v1Tilde,v2Tilde)
            
            if (exist('discretization','var') && strcmp(discretization,'forward'))
                error('Not working')
                grad = generateForwardGradientND( dims,ones(numel(dims),1) );
                
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
            elseif (exist('discretization','var') && strcmp(discretization,'backward'))
                error('Not working')
                grad = generateBackwardGradientND( dims,ones(numel(dims),1) );
                
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
			elseif (exist('discretization','var') && strcmp(discretization,'regularCentral'))
                gradientXf = gradientOperator(size(image1),1,'discretization','forward');
                gradientXb = gradientOperator(size(image1),1,'discretization','backward');
                gradientY = (gradientXf.matrix + gradientXb.matrix)/2;
                
                gradientYf = gradientOperator(size(image1),2,'discretization','forward');
                gradientYb = gradientOperator(size(image1),2,'discretization','backward');
				gradientX = (gradientYf.matrix + gradientYb.matrix)/2;

                
                I2wx = reshape(gradientX * image2(:),size(image2)); %todo implement warp
                I2wy = reshape(gradientY * image2(:),size(image2)); %todo implement warp
                I2w = image2;  %todo implement warp
                
                I1x = 0;
                I1y = 0;
                I2wxx = 0;
                I2wxy = 0;
                I2wyx = 0;
                I2wyy = 0;
                markerOutOfGrid = 0;
			elseif (exist('discretization','var') && strcmp(discretization,'interpolated'))
                
                methodInter = 'spline';

                [M,N] = size(image1);
                
                [gridX,gridY] = meshgrid(1:N,1:M);

                idxx = gridX + reshape(v1Tilde,[M,N]);
                idyy = gridY + reshape(v2Tilde,[M,N]);
                
                markerOutOfGrid = (idxx>=size(idxx,2)) + (idxx<=1) + (idyy>=size(idyy,1)) + (idyy<=1);
                markerOutOfGrid = markerOutOfGrid > 0;

                idxx = max(1,min(N,idxx));
                idxm = max(1,min(N,idxx-0.5));
                idxp = max(1,min(N,idxx+0.5));

                idyy = max(1,min(M,idyy));
                idym = max(1,min(M,idyy-0.5));
                idyp = max(1,min(M,idyy+0.5));
                
                gridXm = max(1,min(N,gridX-0.5));
                gridXp = max(1,min(N,gridX+0.5));
                
                gridYm = max(1,min(M,gridY-0.5));
                gridYp = max(1,min(M,gridY+0.5));
                
                I1x = interp2(image1,gridXp,gridY,methodInter) - interp2(image1,gridXm,gridY,methodInter);
                I1y = interp2(image1,gridX,gridYp,methodInter) - interp2(image1,gridX,gridYm,methodInter);
        
                I2w = interp2(image2,idxx,idyy,methodInter);
                I2wx = interp2(image2,idxp,idyy,methodInter) - interp2(image2,idxm,idyy,methodInter);
                I2wy = interp2(image2,idxx,idyp,methodInter) - interp2(image2,idxx,idym,methodInter);

                %second order derivatives
                I2wxx = (interp2(image2,idxp,idyy,methodInter) + interp2(image2,idxm,idyy,methodInter) - 2*interp2(image2,idxx,idyy,methodInter)) / 2;
                I2wyy = (interp2(image2,idxx,idyp,methodInter) + interp2(image2,idxx,idym,methodInter) - 2*interp2(image2,idxx,idyy,methodInter)) / 2;
                
                I2wxy = (interp2(image2,idxp,idyp,methodInter) - interp2(image2,idxp,idym,methodInter) - (interp2(image2,idxm,idyp,methodInter) - interp2(image2,idxm,idym,methodInter)) ) / 2;
                I2wyx = I2wxy;
                %I2wxx = interp2(I2wx,gridXp,gridY,methodInter) - interp2(I2wx,gridXm,gridY,methodInter);
                %I2wxy = interp2(I2wx,gridX,gridYp,methodInter) - interp2(I2wx,gridX,gridYm,methodInter);
                %I2wyx = interp2(I2wy,gridXp,gridY,methodInter) - interp2(I2wy,gridXm,gridY,methodInter);
                %I2wyy = interp2(I2wy,gridX,gridYp,methodInter) - interp2(I2wy,gridX,gridYm,methodInter);
            else
                error('Not working')
                grad = generateCentralGradientND( dims,ones(numel(dims),1) );
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
            end
            
            
        end
    end
end