classdef nuclearProxDual < basicProx
    properties
    end

    methods
        function obj = nuclearProxDual()
        end

        function applyProx(obj,main,dualNumbers,~)
            for i=1:obj.numVars
                %reshape dual variable back to matrix form
                mat = reshape(main.yTilde{dualNumbers(i)}, main.dims{main.DcP{dualNumbers(i)}});
                
                %retrieve singular values
                [U,S,V] = svd(mat);
                svs = diag(S);
                
                %evaluate L1 prox for singular values
                proxSvs = max(-obj.factor, min(obj.factor, svs));
                
                %reconstruct
                mat = U*diag(proxSvs)*V';
                main.y{dualNumbers(i)} = mat(:);               
            end
        end
    end
end
