classdef flexbox < handle
    
    properties
        primals
        duals
        params
        x
        xOld
        xTilde
        xBar
        
        y
        yOld
        yTilde
        dims
        numberPvars
        PcP %which primal variable belongs to which primal part
        DcP
        DcD
        
        xError
        yError
        
        adaptive
    end
    
    methods
        function obj = flexbox
            % Constructor
            obj.params.tol = 1e-5; 
            obj.params.maxIt = 10000;   
            %obj.params.sigma = 0.1;
            %obj.params.tau = 0.1;
            obj.params.checkError = 100;
            obj.params.theta = 1;
            
            obj.params.showPrimals = 0;
            
            %adaptiveStepsize
            obj.params.adaptivity = 0.5;
            obj.params.delta = 1.5;
            obj.params.eta = 0.95;
            obj.params.s = 255;
            
            %try to use CPP
            obj.params.tryCPP = 0;
            
            obj.primals = {};
            obj.duals = {};
            
            obj.x = {};
            obj.xOld = {};
            obj.xTilde = {};
            obj.xBar = {};
            obj.y = {};
            obj.yOld = {};
            
            obj.dims = {};
            obj.numberPvars = 0;
            
            obj.PcP = {};
            obj.DcP = {};
            obj.DcD = {};
            
            obj.xError = {};
            obj.yError = {};
        end
        
        function number = addPrimalVar(obj,dims)
            %number = addPrimalVar(dims)
            %adds a primal var of dimensions #dims to FlexBox and returns
            %the internal #number
            numberElements = prod(dims);
            
            obj.x{end+1} = zeros(numberElements,1);
            obj.xOld{end+1} = zeros(numberElements,1);
            obj.xBar{end+1} = zeros(numberElements,1);
            obj.xTilde{end+1} = zeros(numberElements,1);
            obj.dims{end+1} = dims;
            obj.numberPvars = obj.numberPvars + 1;
            
            number = obj.numberPvars;
        end
        
        function numberDuals = addTerm(obj,term,corresponding)
            %numberDuals = addTerm(obj,term,corresponding)
            %adds a functional term to FlexBox. The variable #corresponding
            %specifies the internal number of corresponding primal variables.
            %The output #numberDuals contains internal number(s) of the created
            %dual variables or 0 if no dual variable has been created
            numberDuals = 0;

            if (numel(corresponding) ~= term.numPrimals)
                error(['Number of corresponding primals wrong. Expecting ',num2str(term.numPrimals),' variables'])
            end

            %overwrite class name:
            s = superclasses(class(term));

            if (sum(ismember(s, 'primalPart')) > 0)
                %this is a primal part
                obj.primals{end+1} = term;
                obj.PcP{end+1} = corresponding;
            elseif (sum(ismember(s, 'dualPart')) > 0)
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
        end
        
        function showPrimal(obj,number)
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
        
        function doCPP(obj)
            %create function call
            
            mexCallString = [];
            
            %generatePrimalVars
            for i=1:numel(obj.dims)
                mexCallString = [mexCallString,'''primalVar''',',obj.dims{',num2str(i),'},'];
            end
            
            for i=1:numel(obj.primals)
                
                ClassName = class(obj.primals{i});
                
                mexCallString = [mexCallString,'''primal'',','''', ClassName  ,''',','obj.primals{',num2str(i),'}.factor,','obj.PcP{',num2str(i),'},'];
                
                if (strcmp(ClassName,'L1dataTerm') || strcmp(ClassName,'L2dataTerm'))
                    mexCallString = [mexCallString,'obj.primals{',num2str(i),'}.f',','];
                elseif (strcmp(ClassName,'emptyDataTerm'))
                    
                end
            end
            
            for i=1:numel(obj.duals)
                
                ClassName = class(obj.duals{i});
                
                %overwrite class name:
                s = superclasses(ClassName);
                
                if (sum(ismember(s, 'L1IsoProxDual')) > 0 && sum(ismember(s, 'tildeMultiOperatorMultiDual')))
                    ClassName = 'L1dualizedOperatorIso';
                elseif (sum(ismember(s, 'L1AnisoProxDual')) > 0 && sum(ismember(s, 'tildeMultiOperatorMultiDual')))
                    ClassName = 'L1dualizedOperatorAniso';
                elseif (sum(ismember(s, 'L2proxDual')) > 0 && sum(ismember(s, 'tildeMultiOperatorMultiDual')))
                    ClassName = 'L2dualizedOperator';
                end
                
                mexCallString = [mexCallString,'''dual''',',''', ClassName  ,''',','obj.duals{',num2str(i),'}.factor',',','obj.duals{',num2str(i),'}.operator,','obj.DcP{',num2str(i),'},'];
                
                if (strcmp(ClassName,'L1dualizedDataTerm'))
                    mexCallString = [mexCallString,'obj.duals{',num2str(i),'}.f',','];
                end
                %obj.duals{i}
            end
            
            mexCallString = ['flexboxCPP(',mexCallString,'''end''',');']
            %result = flexboxCPP('primal','L1dataTerm',obj.primals{1}.factor,obj.primals{1}.dims,obj.primals{1}.f,'dual','L1gradientAniso',obj.duals{1}.factor,obj.duals{1}.operator,'end');

            
            [resultCPP] = eval(mexCallString);
            figure(1);imagesc(resultCPP,[0,1]);colormap(gray)
        end
        
        function runAlgorithm(obj,varargin)
            %runAlgorithm
            %executes FlexBox and resets the internal iteration counter.
            %The execution can be terminated at any time without loosing
            %the current state
            vararginParser;
            
            if (exist('noInit','var') && noInit == 1)
                
            else
                obj.init(); %init stuff
            end
            
            if (obj.checkCPP())
                disp('C++ support for functional detected');
                disp('Shifting everything to C++ backend');
                obj.doCPP();
            else
                disp('Running flexBox in MATLAB mode');
                reverseStr = [];
                iteration = 1;error = Inf;
                while error>obj.params.tol && iteration < obj.params.maxIt
                    obj.doIteration;

                    if (mod(iteration,obj.params.checkError) == 0)
                        error = obj.calculateError;
                        reverseStr = printToCmd( reverseStr,sprintf(['Iteration: #%d : Residual %.',num2str(-log10(obj.params.tol)),'f'],iteration,error) );
                    end

                    if (obj.params.showPrimals > 0 && mod(iteration,obj.params.showPrimals) == 1)
                        obj.showPrimals
                    end

                    %if (mod(iteration,1) == 0)
                    %    obj.adaptStepsize;
                    %end

                    iteration = iteration + 1;
                end
                printToCmd( reverseStr,'');
            end
        end
    end
    
    methods (Access=protected,Hidden=true )
        %protected methods that can only be accessed from class or
        %subclasses. These methods are hidden!
        
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
            
            %calc xTilde
            for i=1:numel(obj.x)
                obj.xTilde{i} = obj.x{i};
            end
            
            %since xTilde consists of applications of K' to y, it should be
            %done by dual variables
            for i=1:numel(obj.duals)
                obj.duals{i}.xTilde(obj,obj.DcD{i},obj.DcP{i});
            end
            
            %primal prox
            for i=1:numel(obj.primals)
                obj.primals{i}.applyProx(obj,obj.PcP{i});
            end

            %do overrelaxation
            for i=1:numel(obj.x)
                obj.xBar{i} = obj.x{i} + obj.params.theta*(obj.x{i} - obj.xOld{i});
            end
        end
        
        function adaptStepsize(obj)
            [~,p,d] = obj.calculateError;
            
            %if primal residual is massively larger than dual
            if ( p > obj.params.s*d*obj.params.delta )
                for i=1:numel(obj.x)
                    obj.params.tau{i} = obj.params.tau{i} / (1-obj.params.adaptivity);%increase primal steplength
                end
                for i=1:numel(obj.y)
                    obj.params.sigma{i} = obj.params.sigma{i} * (1-obj.params.adaptivity);%decrease dual steplength
                end
                obj.params.adaptivity = obj.params.adaptivity * obj.params.eta;%decrease level of adaptivity
            %if dual residual is massively larger than primal
            elseif (p < obj.params.s*d/obj.params.delta)
                for i=1:numel(obj.x)
                    obj.params.tau{i} = obj.params.tau{i} * (1-obj.params.adaptivity);%decrease primal steplength
                end
                for i=1:numel(obj.y)
                    obj.params.sigma{i} = obj.params.sigma{i} / (1-obj.params.adaptivity);%increase dual steplength
                end
                obj.params.adaptivity = obj.params.adaptivity * obj.params.eta;%decrease level of adaptivity
            end

%             p
%             d
%             obj.params.tau
%             obj.params.sigma
        end
        
        function init(obj)
            %check for primal variables that do not correspond to any term
            %and create an empty term for them
            for i=1:numel(obj.x)
                check = 0;
                for j=1:numel(obj.PcP)
                    if (sum(ismember(obj.PcP{j}, i)) > 0)
                        check = 1;
                    end
                end
                if (~check)
                    obj.addTerm(emptyDataTerm(),i);
                end
            end
            
            %init with zero
            for i=1:numel(obj.x)
                obj.params.tau{i} = 0;
            end
            
            for i=1:numel(obj.y)
                obj.params.sigma{i} = 0;
            end
            
            %init primals
            for i=1:numel(obj.primals)
                obj.primals{i}.init(i,obj);
            end
            
            %init duals
            for i=1:numel(obj.duals)
                obj.duals{i}.init(i,obj);
                
                for j=1:numel(obj.DcP{i})
                    indexTmp = obj.DcP{i}(j);
                    obj.params.tau{ indexTmp } = obj.params.tau{ indexTmp } + obj.duals{i}.myTau{j};
                end
                
                for j=1:numel(obj.DcD{i})
                    indexTmp = obj.DcD{i}(j);
                    obj.params.sigma{ indexTmp } = obj.params.sigma{ indexTmp } + obj.duals{i}.mySigma{j};
                end
            end
            
            %init with zero
            for i=1:numel(obj.x)
                obj.params.tau{i} = 1 ./ obj.params.tau{i};
            end
            
            for i=1:numel(obj.y)
                obj.params.sigma{i} = 1 ./ obj.params.sigma{i};
            end
        end
        
        function [res,resP,resD] = calculateError(obj)
            %calculate first part
            for i=1:numel(obj.x)
                obj.xError{i} = (obj.x{i} - obj.xOld{i}) / obj.params.tau{i};
            end
            for i=1:numel(obj.y)
                obj.yError{i} =  (obj.y{i} - obj.yOld{i}) / obj.params.sigma{i};
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
            if (~obj.params.tryCPP)
                CPPsupport = 0;
            else
                disp('Checking classes for C++ support');
                
                CPPsupport = 1;
                for i=1:numel(obj.primals)
                    if (~obj.primals{i}.CPPsupport)
                        ClassName = class(obj.primals{i});
                        disp(['No support for class ',ClassName]);
                        CPPsupport = 0;
                    end
                end

                for i=1:numel(obj.duals)
                    if (~obj.duals{i}.CPPsupport)
                        ClassName = class(obj.duals{i});
                        disp(['No support for class ',ClassName]);
                        CPPsupport = 0;
                    end
                end
            end
            
            result = CPPsupport;
        end
    end
    
end
