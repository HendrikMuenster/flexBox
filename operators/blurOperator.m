%represents a blur operator
%note: was called convolutionalOperator in an earlier version
classdef blurOperator < basicOperator
    properties
        kernel          %2d kernel
        fftFilter       %kernel in fourier space
        fftFilterC      %complex conjuage of kernel in fourier space
        inputDimension  %dimensions of input
        transposed      %is-transposed-flag in order to decide which kernel version to use
    end

    methods
        function obj = blurOperator(inputDimension,hsize,sigma,varargin)
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
            obj.isMinus = 0;
        end

        function result = mtimes(obj,vector)
            if (obj.transposed)
                result = ifftn(fftn( reshape(vector,obj.inputDimension)  ) .* obj.fftFilterC);
                result = result(:);
            else
                result = ifftn(fftn( reshape(vector,obj.inputDimension)  ) .* obj.fftFilter);
                result = result(:);
            end
            
            if (obj.isMinus)
                result = -result;
            end
        end

        %placeholder
        function result = abs(obj)
            result = 1;
        end

        %placeholder
        function mat = returnMatrix(obj)
            error('Cannot return matrix for convolution operator');
            mat = 1;
        end

        function res = ctranspose(obj)
            res = obj;
            res.transposed = ~res.transposed;
        end

        function result = size(obj,varargin)
            result = prod(obj.inputDimension);

            if (nargin < 2)
                result = [result,result];
            end
        end
        
        function result = uminus(obj)
            result = obj;
            result.isMinus = ~result.isMinus;
        end
        
        function result = getRowSumAbs(obj)
            if (obj.transposed == 1)
                result = sum(abs(obj.kernel(:)));
            else
                result = sum(abs(obj.kernel(:)));
            end
        end

        function result = getMaxRowSumAbs(obj)
            result = max(obj.getRowSumAbs());
        end
    end

end
