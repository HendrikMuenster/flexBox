for i=1:numel(varargin)
    if (mod(i,2)==1)
        assignin('caller',varargin{i}, varargin{i+1})
        if (strcmp(varargin{i},'y1'))
            y1 = varargin{i+1};
        end
        if (strcmp(varargin{i},'y2'))
            y2 = varargin{i+1};
        end
        if (strcmp(varargin{i},'y3'))
            y3 = varargin{i+1};
        end
        if (strcmp(varargin{i},'y4'))
            y4 = varargin{i+1};
        end
        if (strcmp(varargin{i},'v'))
            v = varargin{i+1};
        end
        if (strcmp(varargin{i},'u'))
            u = varargin{i+1};
        end
        if (strcmp(varargin{i},'x'))
            x = varargin{i+1};
        end
        if (strcmp(varargin{i},'Kx'))
            Kx = varargin{i+1};
        end
        if (strcmp(varargin{i},'y'))
            y = varargin{i+1};
        end
        if (strcmp(varargin{i},'Kty'))
            Kty = varargin{i+1};
        end
        if (strcmp(varargin{i},'maxIt'))
            maxIt = varargin{i+1};
        end
        if (strcmp(varargin{i},'displayEvery'))
            displayEvery = varargin{i+1};
        end
        if (strcmp(varargin{i},'tolU'))
            tolU = varargin{i+1};
        end
        if (strcmp(varargin{i},'tolV'))
            tolV = varargin{i+1};
        end
        if (strcmp(varargin{i},'tol'))
            tol = varargin{i+1};
        end
        if (strcmp(varargin{i},'doIsoshrink'))
            doIsoshrink = varargin{i+1};
        end
        if (strcmp(varargin{i},'maxInnerIt'))
            maxInnerIt = varargin{i+1};
        end
        if (strcmp(varargin{i},'maxItOuter'))
            maxItOuter = varargin{i+1};
        end
        if (strcmp(varargin{i},'gtFlow'))
            gtFlow = varargin{i+1};
        end
        if (strcmp(varargin{i},'numsteps'))
            numsteps = varargin{i+1};
        end
        if (strcmp(varargin{i},'steplength'))
            steplength = varargin{i+1};
        end
        if (strcmp(varargin{i},'alphaScaling'))
            alphaScaling = varargin{i+1};
        end
        if (strcmp(varargin{i},'typeNorm'))
            %1 aniso
            %2 iso
            %3 mixed
            typeNorm = varargin{i+1};
        end
        if (strcmp(varargin{i},'stepsize'))
            stepsize = varargin{i+1};
        end
        if (strcmp(varargin{i},'discretization'))
            discretization = varargin{i+1};
        end
        if (strcmp(varargin{i},'useCPP'))
            useCPP = varargin{i+1};
        end
        if (strcmp(varargin{i},'alpha1'))
            alpha1 = varargin{i+1};
        end
        if (strcmp(varargin{i},'numPrimalVars'))
            numPrimalVars = varargin{i+1};
        end
        if (strcmp(varargin{i},'numBreg'))
            numBreg = varargin{i+1};
        end
        if (strcmp(varargin{i},'regularizer'))
            regularizer = varargin{i+1};
        end
        if (strcmp(varargin{i},'dataterm'))
            dataterm = varargin{i+1};
        end
        if (strcmp(varargin{i},'operatorK'))
            operatorK = varargin{i+1};
        end
    end
end

if (~exist('displayEvery','var'))
    displayEvery = 0;
end

if (~exist('doIsoshrink','var'))
    doIsoshrink = 0;
end

