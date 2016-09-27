clear all;close all;clc;

cd 'source';

%%
if(ispc)
    eval(['mex -largeArrayDims -v COMPFLAGS="/openmp $COMPFLAGS" CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp /openmp" LDFLAGS="\$LDFLAGS -fopenmp /openmp" LINKFALGS="$LINKFALGS /openmp" -output ../flexBoxCPP functions.cpp ']);
else
    eval(['mex -largeArrayDims -v  CXXOPTIMFLAGS="-O3 -DNDEBUG" COPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS=" -std=c++11 -fopenmp -O3 -DNDEBUG $CXXFLAGS" LDFLAGS="$LDFLAGS -fopenmp " LINKFALGS="$LINKFALGS -O3 -DNDEBUG " -output ../flexBoxCPP functions.cpp ']);
end
cd '..'