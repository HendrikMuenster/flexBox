%prox for G = alpha / 2 |u_t + \nabla u\cdot v|^2, where v is the unknown
classdef basicOpticalFlow < basicDualizedDataterm
    
    properties
        uxut
        uyut
        uxuy
        ux
        uy
        ut
        ux2
        uy2
        normSquared
        nablaNorm
        breg
    end

    methods
        function obj = basicOpticalFlow(alpha,image1,image2,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            dims = size(image1);
            nPx = prod(dims);
            %obj = obj@primalPart(alpha);%two primal variables
			%obj.numPrimals = 2;
            
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

                idx = repmat([1:N], M,1);
                idy = repmat([1:M]',1,N);

                idxx = idx;
                idyy = idy;
                m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);

                idxx = max(1,min(N,idxx));
                idxm = max(1,min(N,idxx-0.5));
                idxp = max(1,min(N,idxx+0.5));

                idyy = max(1,min(M,idyy));
                idym = max(1,min(M,idyy-0.5));
                idyp = max(1,min(M,idyy+0.5));

                uxStatic = interp2(image2,idxp,idyy,'cubic') - interp2(image2,idxm,idyy,'cubic');
                uxStatic = uxStatic(:);
                uyStatic = interp2(image2,idxx,idyp,'cubic') - interp2(image2,idxx,idym,'cubic');
                uyStatic = uyStatic(:);

                uxStatic(m) = 0.0;
                uyStatic(m) = 0.0;
    
                grad = [uyStatic(:);uxStatic(:)];
                ut = image2(:)-image1(:);
                ut(m) = 0;
            else
                grad = generateCentralGradientND( dims,ones(numel(dims),1) );
                grad = grad*image2(:);
                
                ut = image2(:)-image1(:);
            end
            
            %obj.uy = grad(1:nPx);
            %obj.ux = grad(nPx+1:2*nPx);
            
            A{1} = spdiags(grad(nPx+1:2*nPx),0,nPx,nPx);
            A{2} = spdiags(grad(1:nPx),0,nPx,nPx);
            
            obj = obj@basicDualizedDataterm(alpha,A,-ut(:),varargin);
            
            
            %obj.initStuff;
        end
        
        function initStuff(obj)
            obj.normSquared = max(1e-8,obj.ux.^2+obj.uy.^2);
            
            obj.nablaNorm{1} = obj.ux ./ obj.normSquared;
            obj.nablaNorm{2} = obj.uy ./ obj.normSquared;
            
            obj.uxut = obj.ux.*obj.ut;
            obj.uyut = obj.uy.*obj.ut;
            obj.uxuy = obj.ux .* obj.uy;
            
            obj.ux2 = obj.ux.^2;
            obj.uy2 = obj.uy.^2;
            
            obj.breg{1} = 0;
            obj.breg{2} = 0;
        end

    end
end