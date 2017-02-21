%represents the former base class for dualized terms
%newer classes should directly derive basicTerm
%all data terms that can be dualized should derive this class directly or indirectly
%a dualized data term is of the form \alpha |Au-f|,
%where A is the linear operator, u is the value to be minimized, and f is the data part
%the class initializes required variables for the fixed-point algorithm
classdef basicDualizedDataterm < basicTerm
    methods
        function obj = basicDualizedDataterm(alpha,numPrimals,A,f,varargin)
            if (nargin > 4 && numel(varargin) == 1)
                varargin = varargin{1};
            end
            
            obj = obj@basicTerm(alpha,numPrimals,A,f,varargin);
        end
    end
end
