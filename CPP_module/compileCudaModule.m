cd 'source';

%% Windows Settings:

if ispc
    %this is for visual studio 2013 and must be adjusted for other compilers!
    
    compilerDirectory = 'C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\';
    cudaDirectory = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\';
    
    combilerBin = [compilerDirectory,'bin'];
    compilerIncude = [compilerDirectory,'include'];
    
    CUDA_INCLUDE = [cudaDirectory,'include'];
    CUDA_LIB64 = [cudaDirectory,'\lib\x64'];
end

%% Linux Settings:
if ismac || isunix
    CUDA_INCLUDE='/usr/local/cuda/include';
    CUDA_LIB64='/usr/local/cuda/lib64';
    GCC_LINE = '';
%you may have to change the C++ compiler by adding the following line
%GCC_LINE = 'GCC="/usr/bin/gcc-4.7"';
end

%%

if ismac
    mystring = ['nvcc  -std=c++11 -O3 -m64 -use_fast_math --cuda -arch=compute_50 -code=sm_50 -I"',matlabroot,'/extern/include" functions.cu --output-file "functionsCuda.cpp"'];
elseif isunix
    mystring = ['nvcc  -std=c++11 -O3 -m64 -use_fast_math --cuda -arch=compute_50 -code=sm_50 -I"',matlabroot,'/extern/include" functions.cu --output-file "functionsCuda.cpp"'];
elseif ispc
    mystring = ['nvcc.exe -use_fast_math --cuda -arch=compute_50 -code=sm_50 -ccbin "',combilerBin,'" -I"',matlabroot,'/extern/include" -I"',CUDA_INCLUDE,'" --output-file "functionsCuda.cpp"  "functions.cu"'];
else
    disp('Cannot recognize platform')
end

[status,cmdout] = system(mystring);

dlmwrite('cudaCompileResult.txt',cmdout,'delimiter','')

%%

if ismac
    mex('-largeArrayDims CXXFLAGS="$CXXFLAGS -fopenmp -std=c++11" LDFLAGS="$LDFLAGS -fopenmp" -I"',CUDA_INCLUDE,'" -L"',CUDA_LIB64,'" -L/lib -lcusparse -lcudart -lcufft -lrt -output ../flexBoxCPP -lstdc++ "functionsCuda.cpp"');
elseif isunix
    eval(['mex -v ',GCC_LINE,' -largeArrayDims CXXFLAGS="$CXXFLAGS -fopenmp -fPIC -std=c++11" LDFLAGS="$LDFLAGS -fopenmp" -I"',CUDA_INCLUDE,'" -L"',CUDA_LIB64,'" -L/lib -lcusparse -lcudart -lcufft -lrt -output ../flexBoxCPP -lstdc++ "functionsCuda.cpp"']);
elseif ispc
    eval(['mex -largeArrayDims CXXFLAGS="$CXXFLAGS -fopenmp -std=c++11" LDFLAGS="$LDFLAGS -fopenmp" -I"',CUDA_INCLUDE,'" -L"',CUDA_LIB64,'" -L/lib -lcusparse -lcudart -lcufft -output ../flexBoxCPP "functionsCuda.cpp"']);
end 

delete('functionsCuda.cpp');

cd '..';