function [ reverseStr ] = printToCmd( reverseStr,message )
    message = strrep(message,'%','%%');
    message = strrep(message,'\','\\');

    fprintf([num2str(reverseStr), message]);
    
    %find extra chars
    %nof_extra = length(strfind(message,'%%'));
    %nof_extra = nof_extra + length(strfind(message,'\\'));
    %nof_extra = nof_extra + length(strfind(message,'\n'));

    reverseStr = repmat(sprintf('\b'), 1, length(message));
    
    %% server
%     fid = fopen('log.txt','a');
%     fprintf(fid,message,'\n');
%     fprintf(fid,'\n');
%     fclose(fid);
end

