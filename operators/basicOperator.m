%basic abstract operator class
%deriving classes should override the following functions:
%mtimes: matrix-matrix multiplication, overloaded function, see MATLAB documentation
%returnMatrix: matrix representation of linear operator
%size: size of operator e.g. size of corresponding matrix representation of linear operator
%getMaxRowSumAbs: maximum absolute row sum of matrix representation
classdef basicOperator
    properties
        isMinus
    end

    methods (Abstract)
        abs(obj)
        ctranspose(obj)
        getMaxRowSumAbs(obj)
        getRowSumAbs(obj)
        mtimes(obj,vector)
        returnMatrix(obj)
        size(obj,dim)
    end

    methods
        function result = uminus(obj)
            result = obj;
            result.isMinus = ~result.isMinus;
        end

        function uplus(obj)
            error('You cannot use plus with operator classes. Please use concatOperator instead!')
        end
        

        function plus(obj,B)
            error('You cannot use plus with operator classes. Please use concatOperator instead!')
        end

        function minus(obj,B)
            error('You cannot use minus with operator classes. Please use concatOperator instead!')
        end
    end

end
