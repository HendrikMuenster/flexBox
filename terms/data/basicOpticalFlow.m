% basis class for optical flow data terms u_t + \nabla u\cdot v
classdef basicOpticalFlow < basicDualizedDataterm
    properties
        thisImage1;
        thisImage2;
        termType %brightnessConstancy (default) or gradientConstancy
        constancyDimension
        numSpatial
    end
    
    methods
        function obj = basicOpticalFlow(alpha,image1,image2,varargin)
            if (nargin > 3 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            initVar('discretization','interpolated');
            initVar('vTilde',zeros([size(image1),ndims(image1)]));
            initVar('termType','brightnessConstancy');
            
            if (strcmp(termType,'gradientConstancy'))
                if (~exist('constancyDimension','var'))
                    error('Please specify constancy dimension by adding ''constancyDimension'',''dim''');
                end
            else
                constancyDimension = 0;
            end
            
            
            du = basicOpticalFlow.generateDatatermParts(discretization,image1,image2,vTilde,termType,constancyDimension);
            
            for i=1:numel(du)-1
                A{i} = diagonalOperator(du{i});
            end
            
            obj = obj@basicDualizedDataterm(alpha,numel(A),A,-du{end}(:),varargin);
            
            obj.thisImage1 = image1;
            obj.thisImage2 = image2;
            obj.termType = termType;
            obj.constancyDimension = constancyDimension;
            obj.numSpatial = ndims(image1);
        end
        
        function warpDataterm(obj,vTilde,varargin)
            initVar('discretization','interpolated');
            
            du = obj.generateDatatermParts(discretization,obj.thisImage1,obj.thisImage2,vTilde,obj.termType,obj.constancyDimension);
            
            for i=1:numel(du)-1
                obj.operator{i} = diagonalOperator(du{i});
                obj.operatorT{i} = obj.operator{i};
            end
            
            obj.f{1} = -du{end};
        end
        
    end
    methods(Static)
        function du = generateDatatermParts(discretization,image1,image2,vTilde,termType,constancyDimension)
            %[I1x,I1y,I2w,I2wx,I2wy,I2wxx,I2wxy,I2wyx,I2wyy,markerOutOfGrid] = basicOpticalFlow.generateGradients(discretization,image1,image2,vTilde);
            [I2w,dI2w,markerOutOfGrid] = basicOpticalFlow.generateGradients(discretization,image1,image2,vTilde);
            
            if (strcmp(termType,'brightnessConstancy'))
                %image derivatives are the first du parts
                du = dI2w;
                du{end+1} = I2w - image1;
                
                idx = repmat({':'},1,ndims(vTilde) - 1);
                idx{end+1} = 0;
                for i=1:numel(dI2w)
                    idx{end} = i;
                    
                    du{end} = du{end} - vTilde(idx{:}) .* du{i};
                end
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
                for i=1:numel(du)
                    du{i}(markerOutOfGrid) = 0;
                end
            end
        end
        
        %function [I1x,I1y,I2w,I2wx,I2wy,I2wxx,I2wxy,I2wyx,I2wyy,markerOutOfGrid] = generateGradients(discretization,image1,image2,vTilde)
        function [I2w,dI2w,markerOutOfGrid] = generateGradients(discretization,image1,image2,vTilde)
            
            numSpatial = ndims(image1);
            
			if (exist('discretization','var') && strcmp(discretization,'regularCentral'))
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
                
                for i=1:numSpatial
                    idx{i} = 1:size(image1,i);
                end
                
                grid = cell(1,numSpatial);
                [grid{:}] = ndgrid(idx{:});
                
                %start with empty marker grid
                markerOutOfGrid = zeros(size(grid{1}));
                
                idxI2w = 'image2';
                
                %create grids
                for i=1:numSpatial
                    sizDim = size(image1,i);
                    
                    idx = repmat({':'},1,numSpatial);idx{end+1} = i;
                    
                    gridShift{i} = grid{i} + vTilde(idx{:});
                    gridShiftm{i} = gridShift{i} - 0.5;
                    gridShiftp{i} = gridShift{i} + 0.5;
                    
                    gridm{i} = max(1,min(sizDim,grid{i} - 0.5));
                    gridp{i} = max(1,min(sizDim,grid{i} + 0.5));
                end
                    
                for i=1:numSpatial
                    sizDim = size(image1,i);
                    %create idx for derivative of I2w
                    idx1 = 'image2';
                    idx2 = 'image2';
                    for j=1:numSpatial
                        if (i == j)
                            idx1 = [idx1,',gridShiftp{',num2str(j),'}'];
                            idx2 = [idx2,',gridShiftm{',num2str(j),'}'];
                        else
                            idx1 = [idx1,',gridShift{',num2str(j),'}'];
                            idx2 = [idx2,',gridShift{',num2str(j),'}'];
                        end
                    end
                    idx1 = [idx1,',''spline'''];
                    idx2 = [idx2,',''spline'''];
                    idxI2w = [idxI2w,',gridShift{',num2str(i),'}'];
                    
                    dI2w{i} = eval(['interpn(',idx1,');']) - eval(['interpn(',idx2,');']);
                    
                    check1 = isnan(dI2w{i});
                    check2 = isinf(dI2w{i});

                    if ((sum(check1(:)) + sum(check2(:))) > 0)
                        error('Isinf + Isnan!!!')
                    end

                    markerOutOfGrid = markerOutOfGrid + (grid{i}>=(sizDim)) + (grid{i}<=1);
                    markerOutOfGrid = markerOutOfGrid + (gridShift{i}>=(sizDim)) + (gridShift{i}<=1);
                    markerOutOfGrid = markerOutOfGrid + (gridShiftp{i}>=(sizDim)) + (gridShiftp{i}<=1);
                    markerOutOfGrid = markerOutOfGrid + (gridShiftm{i}>=(sizDim)) + (gridShiftm{i}<=1);
                    
                end
                markerOutOfGrid = markerOutOfGrid > 0;
                
                I2w = eval(['interpn(',idxI2w,',''spline''',');']);
                
                check1 = isnan(I2w);
                check2 = isinf(I2w);

                if ((sum(check1(:)) + sum(check2(:))) > 0)
                    error('Isinf + Isnan!!!')
                end

%                 
%                 idx{3} = 1;
%                 [M,N] = size(image1);
%                 
%                  [gridX,gridY] = meshgrid(1:N,1:M);
% 
%                 
%                 idxx = gridX + reshape(vTilde(:,:,1),[M,N]);
%                 idyy = gridY + reshape(vTilde(:,:,2),[M,N]);
%                 grid{1}
%                 gridShift{1}
%                 gridX
%                 idxx
%                 vTilde(idx{:})
%                 
% 
%                 
%                 markerOutOfGrid = (idxx>=size(idxx,2)) + (idxx<=1) + (idyy>=size(idyy,1)) + (idyy<=1);
%                 markerOutOfGrid = markerOutOfGrid > 0;
% 
%                 idxx = max(1,min(N,idxx));
%                 idxm = max(1,min(N,idxx-0.5));
%                 idxp = max(1,min(N,idxx+0.5));
% 
%                 idyy = max(1,min(M,idyy));
%                 idym = max(1,min(M,idyy-0.5));
%                 idyp = max(1,min(M,idyy+0.5));
%                 
%                 gridXm = max(1,min(N,gridX-0.5));
%                 gridXp = max(1,min(N,gridX+0.5));
%                 
%                 gridYm = max(1,min(M,gridY-0.5));
%                 gridYp = max(1,min(M,gridY+0.5));
%                 
%                 I1x = interp2(image1,gridXp,gridY,methodInter) - interp2(image1,gridXm,gridY,methodInter);
%                 I1y = interp2(image1,gridX,gridYp,methodInter) - interp2(image1,gridX,gridYm,methodInter);
%         
%                 I2w = interp2(image2,idxx,idyy,methodInter);
%                 I2wx = interp2(image2,idxp,idyy,methodInter) - interp2(image2,idxm,idyy,methodInter);
%                 I2wy = interp2(image2,idxx,idyp,methodInter) - interp2(image2,idxx,idym,methodInter);
%                 
%                 I1t = dI2w{i};
%                 I2t = I2wx;
%                 
%                 figure(1);imagesc(dI2w{i});colorbar;
%                 figure(2);imagesc(I2wx);colorbar;
%                 figure(3);imagesc(I1t-I2t);colorbar;
%                 
%                 
                
%                 I1x = interp2(image1,gridXp,gridY,methodInter) - interp2(image1,gridXm,gridY,methodInter);
%                 I1y = interp2(image1,gridX,gridYp,methodInter) - interp2(image1,gridX,gridYm,methodInter);
%         
%                 I2w = interp2(image2,idxx,idyy,methodInter);
%                 I2wx = interp2(image2,idxp,idyy,methodInter) - interp2(image2,idxm,idyy,methodInter);
%                 I2wy = interp2(image2,idxx,idyp,methodInter) - interp2(image2,idxx,idym,methodInter);
% 
%                 %second order derivatives
%                 I2wxx = (interp2(image2,idxp,idyy,methodInter) + interp2(image2,idxm,idyy,methodInter) - 2*interp2(image2,idxx,idyy,methodInter)) / 2;
%                 I2wyy = (interp2(image2,idxx,idyp,methodInter) + interp2(image2,idxx,idym,methodInter) - 2*interp2(image2,idxx,idyy,methodInter)) / 2;
%                 
%                 I2wxy = (interp2(image2,idxp,idyp,methodInter) - interp2(image2,idxp,idym,methodInter) - (interp2(image2,idxm,idyp,methodInter) - interp2(image2,idxm,idym,methodInter)) ) / 2;
%                 I2wyx = I2wxy;
                %I2wxx = interp2(I2wx,gridXp,gridY,methodInter) - interp2(I2wx,gridXm,gridY,methodInter);
                %I2wxy = interp2(I2wx,gridX,gridYp,methodInter) - interp2(I2wx,gridX,gridYm,methodInter);
                %I2wyx = interp2(I2wy,gridXp,gridY,methodInter) - interp2(I2wy,gridXm,gridY,methodInter);
                %I2wyy = interp2(I2wy,gridX,gridYp,methodInter) - interp2(I2wy,gridX,gridYm,methodInter);
            else
                error('Not working');
            end
        end
    end
end