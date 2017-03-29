%represents an empty matrix
classdef zeroOperator < basicOperator
    properties
        nPx;
    end

    methods
        function obj = zeroOperator(nPx1,varargin)
			if (nargin == 1)
				obj.nPx(1) = nPx1;
				obj.nPx(2) = nPx1;
			else
				obj.nPx(1) = nPx1;
				obj.nPx(2) = varargin{1};
			end

        end

        function result = mtimes(~,~)
            result = 0;
        end

        function result = abs(obj)
            result = returnMatrix(obj);
        end

        function mat = returnMatrix(obj)
            mat = sparse(obj.nPx(1),obj.nPx(2));
        end

        function result = size(obj,varargin)
			if (nargin > 1)
                dim = varargin{1};
                result = obj.nPx(dim);
            else
                result = [obj.nPx(1),obj.nPx(2)];
            end
        end

         function res = ctranspose(obj)
             res = obj;
         end
         
        function result = getRowSumAbs(obj)
            result = 0;
        end
        function result = getMaxRowSumAbs(obj)
            result = 0;
        end
    end

end
