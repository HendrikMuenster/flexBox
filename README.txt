% FlexBox
%
% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.04
% Date: 2016-12-16

% FlexBox is copyright ©2016 by Hendrik Dirks and is distributed under the terms of the GNU General Public License (GPL) version 3 (or later).
%
% If you plan to distribute the software (commercially or not). Please contact me for more information.

%%%%%%%%%%%%%%%%
%%% Citation %%%
%%%%%%%%%%%%%%%%

If you use this toolbox please use the following citation
@Article{dirks2015flexbox,
  Title         = {A Flexible Primal-Dual Toolbox},
  Author        = {Dirks, Hendrik},
  Journal       = {ArXiv e-prints},
  Year          = {2016},
  Month         = mar,
  Keywords      = {Mathematics - Optimization and Control, Computer Science - Computer Vision and Pattern Recognition, Computer Science - Mathematical Software, I.4, G.1.6, G.4},
  Primaryclass  = {math.OC}
}
% A preprint of the article can be found at http://arxiv.org/abs/1603.05835

% Examples can be found in the 'examples' folder

% In case you experience any problems, please create an issue at https://github.com/HendrikMuenster/flexBox/issues

%%%%%%%%%%%%%%%%%%%
%%% Compilation %%%
%%%%%%%%%%%%%%%%%%%

C++:
To compile the C++ module, simply run the MATLAB script "compileMexModule.m" in the folder "CPP_module"

CUDA:
To compile the Cuda module, open the "compileCudaModule.m" and change the given PATH variables according to your computer. Afterwards, run "compileCudaModule.m".

Requirements:
The script assumes, that you have installed a full version of Microsoft Visual Studio 20XX, which is connected to your CUDA SDK. The Visual Studio version must coincide with your Matlab standard C++ compiler. 

Issues for University of Münster users:
There are some known problems for CUDA computers at the University of Münster. You may have to preload libraries. To do so, simply start Matlab with:
env LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/cuda-7.5/lib64/libcudart.so.7.5 /usr/local/cuda-7.5/lib64/libcusparse.so.7.5" matlab

%%%%%%%%%%%%%%%%%
%%% Changelog %%%
%%%%%%%%%%%%%%%%%
16.12.2016 (Version 1.04)
- new type of terms 'constraint' added. For now, box-constraint, non-neg. and non-pos constraints are available
- fused the basicDualizedDataterm.m and basicDualizedOperator.m. An operator is now a data term with f=0
- new MATLAB->CPP wrapper
- several bugfixes

14.07.2016 (Version 1.03)
- CUDA module added

07.04.2016 (Version 1.02)
- Added FFT-based convolution operator, downsampling and concatenation operator
- Added examples for deblurring and operator concatenation
- Speed improvements for C++ version
31.03.2016 (Version 1.01)
- generalized the concept of matrices to linear operators and created classes for gradient, diagonal, identity and zero operator for MATLAB and C++ version to improve speed and flexibility