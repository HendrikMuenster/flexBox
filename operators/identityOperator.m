%represents an identity matrix
classdef identityOperator < basicOperator
    properties
        nPx;
    end

    methods
        function obj = identityOperator(nPx,varargin)
            obj.nPx = nPx;
            obj.isMinus = 0;
        end

        function result = mtimes(obj,vector)
            if (obj.isMinus)
                result = -vector;
            else
                result = vector;
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
        end

        function result = size(obj,varargin)
            if (nargin < 2)
                result = [obj.nPx,obj.nPx];
            else
                result = obj.nPx;
            end
        end

%         function res = ctranspose(obj)
%             res = obj;
%             res.transposed = ~obj.transposed;
%         end

        function result = getMaxRowSumAbs(obj)
            result = 1; %matrix representation is identiy matrix -> max absolute value is 1
        end
    end

end
