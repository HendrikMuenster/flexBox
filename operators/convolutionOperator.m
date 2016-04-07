%represents a convolution operator
classdef convolutionOperator < basicOperator
    properties
        kernel
        fftFilter
        fftFilterC
        inputDimension
        transposed
    end
    
    methods
        function obj = convolutionOperator(inputDimension,hsize,sigma,varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;
            
            if (mod(hsize,2) == 0)
                error('kernel size must be odd')
            end
            
            if (sigma < 0)
                error('sigma size must positive')
            end

            obj.inputDimension = inputDimension;

            obj.kernel = fspecial('gaussian', hsize, sigma);
            
            obj.fftFilter = zeros(obj.inputDimension);
            obj.fftFilter(1:size(obj.kernel,1),1:size(obj.kernel,2)) = obj.kernel;
            
            %center kernel
            obj.fftFilter = circshift(obj.fftFilter,-(size(obj.kernel,1)-1)/2,1);
            obj.fftFilter = circshift(obj.fftFilter,-(size(obj.kernel,2)-1)/2,2);
            
            %precalculate fft
            obj.fftFilter = fftn(obj.fftFilter);
            obj.fftFilterC = conj(obj.fftFilter);
            
            obj.transposed = 0;
        end
        
        function result = mtimes(obj,vector)
            if (obj.transposed)
                result = ifftn(fftn( reshape(vector,obj.inputDimension)  ) .* obj.fftFilter);
                result = result(:);
            else
                result = ifftn(fftn( reshape(vector,obj.inputDimension)  ) .* obj.fftFilterC);
                result = result(:);
            end
            
        end
        
        function result = abs(obj)
            result = 1;
        end
        
        function mat = returnMatrix(obj)
            mat = 1;
        end
        
        function res = ctranspose(obj)
            res = obj;
            res.transposed = ~res.transposed;
        end
        
        function result = size(obj,dim)
            result = prod(obj.inputDimension);
        end
        
        function result = getMaxRowSumAbs(obj)
            if (obj.transposed == 1)
                result = sum(abs(obj.kernel(:)));
            else
                result = sum(abs(obj.kernel(:)));
            end
        end
    end
    
end

