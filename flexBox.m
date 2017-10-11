classdef flexBox < handle

    properties
        params      %FlexBox params (tolerances, number of iterations, etc. see documentation)

        duals       %cell array of dual terms

        x           %current iteration values for primal part
        xOld        %old iteration values for primal part
        xBar        %overrelaxed value for x

        y           %current iteration values for dual part
        yOld        %old iteration values for dual part
        yTilde      %intermediate value for primal dual algortihm

        dims        %dimensionals for primal variables
        numberPvars %internal unique identifiers for primal variables
        DcP         %which primal variable(s) belongs to which dual part
        DcD         %which dual variables(s) belong to which dual part

        xError      %error estimates for primal parts
        yError      %error estimates for dual parts
        error       %combined error estimate
        
        firstRun    %internal flag specifying if algortihm has already been run
    end
    
    methods
        function obj = flexBox
            % Constructor
            obj.params.tol = 1e-5;
            obj.params.maxIt = 10000;
            obj.params.checkError = 100;
            obj.params.theta = 1;
            obj.params.verbose = 0;

            obj.params.showPrimals = 0;

            obj.params.tryCPP = 0;
            obj.params.relativePathToMEX = 'flexBox_CPP/';

            obj.duals = {};

            obj.x = {};
            obj.xOld = {};
            obj.xBar = {};
            obj.y = {};
            obj.yOld = {};

            obj.dims = {};
            obj.numberPvars = 0;

            obj.DcP = {};
            obj.DcD = {};

            obj.xError = {};
            obj.yError = {};
            
            obj.firstRun = true;
        end

        function number = addPrimalVar(obj,dims)
            %number = addPrimalVar(dims)
            %adds a primal var of dimensions #dims to FlexBox and returns
            %the internal #number

            if (isscalar(dims))
                dims = [dims,1];
            end

            numberElements = prod(dims);

            obj.x{end+1} = zeros(numberElements,1);
            obj.xOld{end+1} = zeros(numberElements,1);
            obj.xBar{end+1} = zeros(numberElements,1);
            obj.dims{end+1} = dims;
            obj.numberPvars = obj.numberPvars + 1;

            number = obj.numberPvars;
        end

        function numberDuals = addTerm(obj,term,corresponding)
            %numberDuals = addTerm(obj,term,corresponding)
            %adds a functional term to FlexBox. The variable #corresponding
            %specifies the internal number of corresponding primal variables.
            %The output #numberDuals contains internal number(s) of the created
            %dual variables
            if (numel(corresponding) ~= term.numPrimals)
                error(['Number of corresponding primals wrong. Expecting ',num2str(term.numPrimals),' variables'])
            end

            %this is a dual part
            obj.duals{end+1} = term;
            %save corresponding Dual->Primal,Dual->Dual
            obj.DcP{end+1} = corresponding;
            obj.DcD{end+1} = numel(obj.y) + 1 : numel(obj.y) + term.numVars;

            numberDuals = obj.DcD{end};

            %create variables
            for i=1:term.numVars
                obj.y{end+1} = zeros(term.length{i},1);
                obj.yOld{end+1} = zeros(term.length{i},1);
            end

        end

        function showPrimal(obj,number)
          %helper function to display primal with internal number #number as 2D
          %or 3D image
            if (numel(obj.dims{number}) == 2)
                figure(number);clf;imagesc(reshape(obj.x{number},obj.dims{number}),[0,1]);axis image;colormap(gray)
            elseif (numel(obj.dims{number}) == 3)
                dims2 = obj.dims{number}(1:2);
                pixelPerSlice = prod(obj.dims{number}(1:2));
                showSlice = 1;

                indexStart = (showSlice-1)*pixelPerSlice + 1;
                endStart = (showSlice)*pixelPerSlice;

                figure(number);clf;imagesc(reshape(obj.x{number}(indexStart:endStart),dims2),[0,1]);axis image;colormap(gray)
            end
        end

        function showPrimals(obj)
          %displays all primal variables (if 2D or 3D)
            for i=1:numel(obj.x)
                obj.showPrimal(i);
            end
            drawnow;
        end

        function u = getPrimal(obj,number)
            %u = getPrimal(number)
            %returns the primal variable specified by #number and reshapes the
            %variable to the correct size
            u = reshape(obj.x{number},obj.dims{number});
        end

        function y = getDual(obj,number)
            %y = getDual(number)
            %returns the dual variable specified by #number as a vector
            y = obj.y{number};
        end

        function runAlgorithm(obj,varargin)
            %runAlgorithm
            %executes FlexBox and resets the internal iteration counter.
            %The execution can be terminated at any time without loosing
            %the current state
            vararginParser;

            if (exist('noInit','var') && noInit == 1)

            else
				%initialize tau and sigma
                obj.init();
            end

			%check if C++ module is activated and compiled
            if (obj.checkCPP())
                obj.doCPP();
            else
                reverseStr = [];

                if obj.firstRun
                    obj.error = Inf;
                    obj.firstRun = false;
                else
                    obj.error = obj.calculateError();
                end
                
                iteration = 1;
                while obj.error > obj.params.tol && iteration <= obj.params.maxIt
                    obj.doIteration;

                    if (mod(iteration,obj.params.checkError) == 0)
                        obj.error = obj.calculateError();
                        reverseStr = printToCmd( reverseStr, sprintf(['Iteration: #%d : Residual %.', num2str(-log10(obj.params.tol)), 'f', '\n'], iteration, obj.error) );
                    end

                    if (obj.params.showPrimals > 0 && mod(iteration,obj.params.showPrimals) == 1)
                        obj.showPrimals;
                    end

                    iteration = iteration + 1;
                end
                printToCmd( reverseStr,'');
            end
        end
    end

    methods (Access=protected,Hidden=true )
        %protected methods that can only be accessed from class or
        %subclasses. These methods are hidden!

        function doCPP(obj)
            %create function call
            [resultCPP{1:numel(obj.x)+numel(obj.y)}] = eval('flexBoxCPP(obj);');

            for i=1:numel(obj.x)
                obj.x{i} = resultCPP{i};
            end

            for i=1:numel(obj.y)
                obj.y{i} = resultCPP{numel(obj.x)+i};
            end

            obj.firstRun = false;
        end

            function doIteration(obj)
            %save old
            for i=1:numel(obj.xOld)
                obj.xOld{i} = obj.x{i};
            end

            for i=1:numel(obj.yOld)
                obj.yOld{i} = obj.y{i};
            end

            %calc yTilde and prox
            for i=1:numel(obj.duals)
                %input are numbers of the dual variables they shall update
                %and the primal variables they correspond to
                obj.duals{i}.yTilde(obj,obj.DcD{i},obj.DcP{i});
                obj.duals{i}.applyProx(obj,obj.DcD{i},obj.DcP{i});
            end

            %primal update is x = x - K^ty
            for k=1:numel(obj.duals)
                dualNumbers = obj.DcD{k};
                primalNumbers = obj.DcP{k};
                for i=1:numel(dualNumbers)
                    for j=1:numel(primalNumbers)
                        operatorNumber = numel(primalNumbers)*(i-1) + j;
                        obj.x{primalNumbers(j)} = obj.x{primalNumbers(j)} - obj.params.tau{primalNumbers(j)}.*(obj.duals{k}.operatorT{operatorNumber} * obj.y{dualNumbers(i)});
                    end
                end
            end

            %do overrelaxation
            for i=1:numel(obj.x)
                obj.xBar{i} = obj.x{i} + obj.params.theta*(obj.x{i} - obj.xOld{i});
            end
        end

%         function adaptStepsize(obj)
%             [~,p,d] = obj.calculateError;
%
%             %if primal residual is massively larger than dual
%             if ( p > obj.params.s*d*obj.params.delta )
%                 for i=1:numel(obj.x)
%                     obj.params.tau{i} = obj.params.tau{i} / (1-obj.params.adaptivity);%increase primal steplength
%                 end
%                 for i=1:numel(obj.y)
%                     obj.params.sigma{i} = obj.params.sigma{i} * (1-obj.params.adaptivity);%decrease dual steplength
%                 end
%                 obj.params.adaptivity = obj.params.adaptivity * obj.params.eta;%decrease level of adaptivity
%             %if dual residual is massively larger than primal
%             elseif (p < obj.params.s*d/obj.params.delta)
%                 for i=1:numel(obj.x)
%                     obj.params.tau{i} = obj.params.tau{i} * (1-obj.params.adaptivity);%decrease primal steplength
%                 end
%                 for i=1:numel(obj.y)
%                     obj.params.sigma{i} = obj.params.sigma{i} / (1-obj.params.adaptivity);%increase dual steplength
%                 end
%                 obj.params.adaptivity = obj.params.adaptivity * obj.params.eta;%decrease level of adaptivity
%             end
%
% %             p
% %             d
% %             obj.params.tau
% %             obj.params.sigma
%         end

        function init(obj)
            %init tau and sigma with all zero vectors
            for i=1:numel(obj.x)
                obj.params.tau{i} = zeros(numel(obj.x{i}),1);
            end

            for i=1:numel(obj.y)
                obj.params.sigma{i} = zeros(numel(obj.y{i}),1);
            end

            %init duals
            for i=1:numel(obj.duals)
                obj.duals{i}.init();

                %sum up tau
                for j=1:numel(obj.DcP{i})
                    indexTmp = obj.DcP{i}(j);

                    obj.params.tau{ indexTmp } = obj.params.tau{ indexTmp } + obj.duals{i}.myTau{j};
                end

                %sum up sigma
                for j=1:numel(obj.DcD{i})
                    indexTmp = obj.DcD{i}(j);
                    obj.params.sigma{ indexTmp } = obj.params.sigma{ indexTmp } + obj.duals{i}.mySigma{j};
                end
            end

            %calculate reciprocals
            for i=1:numel(obj.x)
                obj.params.tau{i} = 1 ./ max(0.0001,obj.params.tau{i});
            end

            for i=1:numel(obj.y)
                obj.params.sigma{i} = 1 ./ max(0.0001,obj.params.sigma{i});
            end
        end

        function [res,resP,resD] = calculateError(obj)
            %calculates residual in primal dual algorithm

            %calculate first part
            for i=1:numel(obj.x)
                obj.xError{i} = (obj.x{i} - obj.xOld{i}) ./ obj.params.tau{i};
            end
            for i=1:numel(obj.y)
                obj.yError{i} =  (obj.y{i} - obj.yOld{i}) ./ obj.params.sigma{i};
            end

            %calculate second part
            for i=1:numel(obj.duals)
                obj.duals{i}.xError(obj,obj.DcD{i},obj.DcP{i});
                obj.duals{i}.yError(obj,obj.DcD{i},obj.DcP{i});
            end

            %sum up
            resP = 0;
            resD = 0;
            for i=1:numel(obj.x)
                respTMP = sum(abs(obj.xError{i}));
                resP = resP + respTMP;
                %resPList{i} = respTMP;
            end
            resP = resP / numel(obj.x);

            for i=1:numel(obj.y)
                resDTMP = sum(abs(obj.yError{i}));
                resD = resD + resDTMP;
                %resDList{i} = respTMP;
            end
            resD = resD / numel(obj.y);

            elements = 0;
            for i=1:numel(obj.dims)
                elements = elements + prod(obj.dims{i});
            end

            res = (resD + resP) / elements;
        end

        function result = checkCPP(obj)
            %checkCPP
            %checks if the tryCPP parameter is true and, if yes, checks if the C++ module is compiled
			%if the C++ module is compiled, but cannot be found the folder specified by obj.params.relativePathToMEX, a warning is displayed
            if (~obj.params.tryCPP)
                CPPsupport = 0;
            elseif (obj.params.tryCPP)
                absPathToMEX = strcat(fileparts(mfilename('fullpath')), '/', obj.params.relativePathToMEX);
                if (exist(absPathToMEX, 'dir') ~= 7) %dir is not correct. Try to find it through path
                    CPPsupport = 0;
                    disp(['Warning: relative Path to MEX-File is not correct! The default path is stored in params.relativePathToMEX']);
                else
                    %make sure the intended MEX file is called
                    addpath(absPathToMEX);
                end

                if (exist('flexBoxCPP','file') ~= 3)
                    CPPsupport = 0;
                    disp(['Warning: C++ module is not compiled!']);
                    disp(['Running in MATLAB mode']);
                else
                    CPPsupport = 1;
					if (obj.params.verbose > 0)
						disp(['using MEX-File: ', which('flexBoxCPP')]);
						disp('Running in C++ mode');
					end
                end
            end
            result = CPPsupport;
        end
    end

end
