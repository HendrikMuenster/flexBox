clear all;close all;clc;

cd 'source';

listOfFiles = dir('.');

cppFileString = '';
for i=1:numel(listOfFiles)
    if (~listOfFiles(i).isdir && ~strcmp(listOfFiles(i).name,'functions.cpp'))
        cppFileString = [cppFileString,' ',listOfFiles(i).name];
    end
end


%%
if(ispc)
    eval(['mex -largeArrayDims -v COMPFLAGS="/openmp $COMPFLAGS" CXXFLAGS="\$CXXFLAGS -fopenmp /openmp" CFLAGS="\$CFLAGS -fopenmp /openmp" LDFLAGS="\$LDFLAGS -fopenmp /openmp" LINKFALGS="$LINKFALGS /openmp" -output ../flexBoxCPP functions.cpp ']);
else
    eval(['mex -largeArrayDims -v COMPFLAGS=" $COMPFLAGS" CXXFLAGS="\$CXXFLAGS -fopenmp " CFLAGS="\$CFLAGS -fopenmp " LDFLAGS="\$LDFLAGS -fopenmp " LINKFALGS="$LINKFALGS " -output ../flexBoxCPP functions.cpp ']);

end
cd '..'