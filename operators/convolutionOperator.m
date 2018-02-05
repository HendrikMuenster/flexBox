%represents a convolution operator
classdef convolutionOperator < basicOperator
    properties
        kernel          %2d kernel
        fftFilter       %kernel in fourier space
        fftFilterC      %complex conjuage of kernel in fourier space
        inputDimension  %dimensions of input
        transposed      %is-transposed-flag in order to decide which kernel version to use
    end

    methods
        function obj = convolutionOperator(kernel, inputDimension, varargin)
            if (nargin > 2 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            obj.kernel = kernel; 
            obj.inputDimension = inputDimension;

            %precalculate fft with truncation and padding
            obj.fftFilter = fftn(obj.kernel, inputDimension);
            obj.fftFilterC = conj(obj.fftFilter);

            obj.transposed = 0;
            obj.isMinus = 0;
        end

        function result = mtimes(obj,vector)
            if (obj.transposed)
                result = ifftn(fftn(reshape(vector,obj.inputDimension), obj.inputDimension) .* obj.fftFilterC);
                result = result(:);
            else
                result = ifftn(fftn(reshape(vector,obj.inputDimension), obj.inputDimension) .* obj.fftFilter);
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
