%represents an identity matrix
classdef identityOperator < basicOperator
    properties
        nPx;
    end

    methods
        function obj = identityOperator(nPx1,varargin)
			if (nargin == 1)
				obj.nPx(1) = nPx1;
				obj.nPx(2) = nPx1;
			else
				obj.nPx(1) = nPx1;
				obj.nPx(2) = varargin{1};
			end
            obj.isMinus = 0;
        end

        function result = mtimes(obj,vector)
            if (obj.isMinus)
                result = -vector;
            else
                result = vector;
            end

      			if (obj.nPx(1) < obj.nPx(2))
      				result = result(1:obj.nPx(1));
      			elseif (obj.nPx(1) > obj.nPx(2))
      				result(obj.nPx(2)+1:obj.nPx(1)) = 0;
      			end
        end

        function result = abs(obj)
            if (obj.isMinus)
                result = -returnMatrix(obj);
            else
                result = returnMatrix(obj);
            end

        end

        function result = returnMatrix(obj)
            if (obj.isMinus)
                result = -speye(obj.nPx);
            else
                result = speye(obj.nPx);
            end

			result = result(1:obj.nPx(1),1:obj.nPx(2));
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
            result = 1; %matrix representation is identiy matrix -> max absolute value is 1
         end

        function result = getMaxRowSumAbs(obj)
            result = 1; %matrix representation is identiy matrix -> max absolute value is 1
        end
    end

end
