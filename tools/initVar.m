function initVar(var,default )
	callerVars = evalin('caller','whos');

    check = 0;
    for i=1:numel(callerVars)
        if (strcmp(callerVars(i).name,var))
            check = 1;
        end
    end
    
    if (~check)
        assignin('caller',var,default)
    end
end