%term for alpha / 2 |u_t + \nabla u\cdot v|_2^2, where v is the unknown
classdef L2opticalFlowTerm < basicOpticalFlow & L2DataProxDual
    methods
        function obj = L2opticalFlowTerm(alpha,image1,image2,varargin)
            obj = obj@basicOpticalFlow(alpha,image1,image2,varargin);
        end
    end
end