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
% Version 1.1
% Date: 2017-04-04

% FlexBox is copyright ©2016-2017 by Hendrik Dirks
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

Download the C++ submodule and create compile files using cmake.

%%%%%%%%%%%%%%%%%%%
%%%%%% Issues %%%%%
%%%%%%%%%%%%%%%%%%%
Issues for University of Münster users:
There are some known problems for CUDA computers at the University of Münster. You may have to preload libraries. To do so, simply start Matlab with:
env LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/cuda-7.5/lib64/libcudart.so.7.5 /usr/local/cuda-7.5/lib64/libcusparse.so.7.5" matlab
