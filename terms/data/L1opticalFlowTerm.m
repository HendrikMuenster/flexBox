%term representing alpha |u_t + \nabla u\cdot v|_1, where v is the unknown
classdef L1opticalFlowTerm < basicOpticalFlow & L1DataProxDual

    methods
        function obj = L1opticalFlowTerm(alpha,image1,image2,varargin)
            obj = obj@basicOpticalFlow(alpha,image1,image2,varargin);
        end
    end
end