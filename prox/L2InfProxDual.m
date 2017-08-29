classdef L2InfProxDual < basicProx
    properties
    end

    methods
        function obj = L2InfProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            yTildeNorm = 0;
            for i=1:obj.numVars
                yTildeNorm = yTildeNorm + main.yTilde{dualNumbers(i)}.^2;
            end
            yTildeNorm = sqrt(yTildeNorm);


            h = @(lambda) (yTildeNorm > lambda) .* (yTildeNorm - lambda);
            g = @(lambda) sum(h(lambda)) - obj.factor;
            %h = @(lambda) (sortyTildeNorm > lambda) .* (sortyTildeNorm - lambda);
            %g = @(lambda) sum(h(lambda)) - obj.factor;
            %lambda = max(0,fzero(g, max(yTildeNorm)));

            sortyTildeNorm = sort(yTildeNorm, 'descend');
            for findIndex=1:length(sortyTildeNorm)
                if g(sortyTildeNorm(findIndex)) > 0
                    break
                end
            end
       
            lambda = 0;
            if ~isempty(findIndex) && findIndex ~= 1 && (findIndex ~= length(sortyTildeNorm) || g(sortyTildeNorm(findIndex)) > 0)
                lambda = (sum(sortyTildeNorm(1:findIndex-1)) - obj.factor) / (findIndex - 1);
            end

            %disp(findIndex)
            %disp(lambda)
            %disp(sortyTildeNorm)
            %disp('-------------')
            %waitforbuttonpress;

            %disp(lambda);
            %zeros = linspace(-1,max(yTildeNorm) + 1, 100);
            %if length(zeros) > 3
            %    disp("yeah")
            %end
            %gZeros = g(zeros);
            %h = figure(1);
            %plot(zeros, gZeros);
            %hold on;
            %plot([min(zeros)-10, max(zeros)+10], [0, 0]);
            %waitforbuttonpress;
            %close(h);
            %if sum(yTildeNorm > lambda) > 2
            %    disp("hooray")
            %end




            for i=1:obj.numVars
                main.y{dualNumbers(i)} = (yTildeNorm > lambda) .* main.yTilde{dualNumbers(i)} .* (1 - lambda ./ yTildeNorm);
                main.y{dualNumbers(i)}(yTildeNorm <= lambda) = 0;
            end



%             %first guess for main.y
%             for i=1:obj.numVars
%                 main.y{dualNumbers(i)} = main.yTilde{dualNumbers(i)};
%             end
%             yNorm = yTildeNorm;
%
%             counter = 1;
%             while (sum(yNorm) >= obj.factor) && (counter <= length(sortyTildeNorm))
%                 lambda = sortyTildeNorm(counter);
%
%                 %update guess for main.y
%                 for i=1:obj.numVars
%                     main.y{dualNumbers(i)} = (yTildeNorm > lambda) .* main.yTilde{dualNumbers(i)} .* (1 - lambda ./ yTildeNorm);
%                     main.y{dualNumbers(i)}(yTildeNorm <= lambda) = 0;
%                 end
%
%                 yNorm = 0;
%                 for i=1:obj.numVars
%                     yNorm = yNorm + main.y{dualNumbers(i)}.^2;
%                 end
%                 yNorm = sqrt(yNorm);
%                 %yNorm = yTildeNorm - lambda;
%
%                 counter = counter + 1;
%             end
        end
    end
end
