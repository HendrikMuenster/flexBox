%represents base class for terms containing the operator
%K(u,v) = [u_x + v_y;u_y - v_x]
%corresponds to two primal variables (u,v) and
classdef basicDivCurl < dualPart
    properties
    end

    methods
        function obj = basicDivCurl(alpha,dims,varargin)
            if (nargin > 2 == numel(varargin) == 1)
                varargin = varargin{1};
            end
            vararginParser;

            if (exist('discretization','var') && strcmp(discretization,'backward'))
                opTmp = generateBackwardGradientND( dims,ones(numel(dims),1) );
            else
                opTmp = generateForwardGradientND( dims,ones(numel(dims),1) );
            end

            obj = obj@dualPart(alpha);
            obj.numVars = 2; %divergence produces scalar quantity
            for i=1:numel(dims)
                obj.length{i} = prod(dims);
                obj.operator{i} = opTmp( (i-1)*prod(dims) + 1 : i * prod(dims),: )';
                obj.myTau{i} = 2;
            end
            obj.mySigma{1} = 2*numel(dims);
            obj.mySigma{2} = 2*numel(dims);

            obj.operator{3} = obj.operator{2};
            obj.operator{4} = -obj.operator{1};
        end

        function yTilde(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                main.yTilde{dualNumbers(i)} = main.y{dualNumbers(i)};
            end

            for i=1:numel(dualNumbers)
                for j=1:2
                    operatorNumber = 2*(i-1) + j;
                    main.yTilde{dualNumbers(i)} = main.yTilde{dualNumbers(i)} + main.params.sigma{dualNumbers(i)} * (obj.operator{operatorNumber}* main.xBar{primalNumbers(j)});
                end
            end
        end

        function xTilde(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(primalNumbers)
                for j=1:2
                    operatorNumber = 2*(j-1) + i;
                    main.xTilde{primalNumbers(i)} = main.xTilde{primalNumbers(i)} - main.params.tau{primalNumbers(i)}*(obj.operator{operatorNumber}' * main.y{dualNumbers(j)});
                end
            end
        end

        function yError(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(dualNumbers)
                for j=1:2
                    operatorNumber = 2*(i-1) + j;
                    main.yError{dualNumbers(i)} = main.yError{dualNumbers(i)} - obj.operator{operatorNumber}* (main.xOld{primalNumbers(j)}-main.x{primalNumbers(j)});
                end
            end
        end

        function xError(obj,main,dualNumbers,primalNumbers)
            for i=1:numel(primalNumbers)
                for j=1:2
                    operatorNumber = 2*(j-1) + i;
                    main.xError{primalNumbers(i)} = main.xError{primalNumbers(i)} - obj.operator{operatorNumber}' * (main.yOld{dualNumbers(j)}-main.y{dualNumbers(j)});
                end
            end
        end
    end
end
