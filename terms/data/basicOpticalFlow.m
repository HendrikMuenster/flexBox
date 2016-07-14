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
            if (nargin > 2 && numel(varargin) == 1)
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
            
            obj = obj@basicDualizedDataterm(alpha,A,-ut(:),varargin);
            
            obj.thisImage1 = image1;
            obj.thisImage2 = image2;
            obj.termType = termType;
            obj.constancyDimension = constancyDimension;
        end
        
        function warpDataterm(obj,v1Tilde,v2Tilde,varargin)
            initVar('discretization','interpolated');
            
            [ux,uy,ut] = obj.generateDatatermParts(discretization,obj.thisImage1,obj.thisImage2,v1Tilde,v2Tilde,obj.termType,obj.constancyDimension);
            
            obj.operator{1} = diagonalOperator(ux);
            obj.operatorT{1} = diagonalOperator(ux);
            
            obj.operator{2} = diagonalOperator(uy);
            obj.operatorT{2} = diagonalOperator(uy);
            
            obj.f = -ut;
            
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
            
            ux(markerOutOfGrid) = 0;
            uy(markerOutOfGrid) = 0;
            ut(markerOutOfGrid) = 0;
        end
        
        function [I1x,I1y,I2w,I2wx,I2wy,I2wxx,I2wxy,I2wyx,I2wyy,markerOutOfGrid] = generateGradients(discretization,image1,image2,v1Tilde,v2Tilde)
            
            if (exist('discretization','var') && strcmp(discretization,'forward'))
                grad = generateForwardGradientND( dims,ones(numel(dims),1) );
                
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
            elseif (exist('discretization','var') && strcmp(discretization,'backward'))
                grad = generateBackwardGradientND( dims,ones(numel(dims),1) );
                
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
			elseif (exist('discretization','var') && strcmp(discretization,'interpolated'))

                [M,N] = size(image1);
                
                [gridX,gridY] = meshgrid(1:N,1:M);

                idxx = gridX + v1Tilde;
                idyy = gridY + v2Tilde;
                
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
                
                I1x = interp2(image1,gridXp,gridY,'bicubic') - interp2(image1,gridXm,gridY,'bicubic');
                I1y = interp2(image1,gridX,gridYp,'bicubic') - interp2(image1,gridX,gridYm,'bicubic');
        
                I2w = interp2(image2,idxx,idyy,'bicubic');
                I2wx = interp2(image2,idxp,idyy,'bicubic') - interp2(image2,idxm,idyy,'bicubic');
                I2wy = interp2(image2,idxx,idyp,'bicubic') - interp2(image2,idxx,idym,'bicubic');

                %second order derivatives
                I2wxx = interp2(I2wx,gridXp,gridY,'bicubic') - interp2(I2wx,gridXm,gridY,'bicubic');
                I2wxy = interp2(I2wx,gridX,gridYp,'bicubic') - interp2(I2wx,gridX,gridYm,'bicubic');
                I2wyx = interp2(I2wy,gridXp,gridY,'bicubic') - interp2(I2wy,gridXm,gridY,'bicubic');
                I2wyy = interp2(I2wy,gridX,gridYp,'bicubic') - interp2(I2wy,gridX,gridYm,'bicubic');
            else
                grad = generateCentralGradientND( dims,ones(numel(dims),1) );
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
            end
            
            
        end
    end
end