%basic abstract operator class
%deriving classes should override the following functions:
%mtimes: matrix-matrix multiplication, overloaded function, see MATLAB documentation
%returnMatrix: matrix representation of linear operator
%size: size of operator e.g. size of corresponding matrix representation of linear operator
%getMaxRowSumAbs: maximum absolute row sum of matrix representation
classdef basicOperator
    properties
    end

    methods (Abstract)
        mtimes(obj,vector)
        returnMatrix(obj)
        size(obj,dim)
        getMaxRowSumAbs(obj)
    end

end
