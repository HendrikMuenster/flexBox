# FlexBox - A **Flex**ible Primal-Dual Tool**Box**

## Introduction
**FlexBox** is a flexible MATLAB toolbox for finite dimensional convex variational problems in image processing and beyond

Nowadays, many problems in image processing consist of minimizing convex energies. Typically these problems can be written as

![Problem formulation][generalPrimalFormulation] ,

where ![A](https://latex.codecogs.com/svg.latex?A) denotes a linear operator and
![F](https://latex.codecogs.com/svg.latex?F) is a proper, convex and lower-semicontinuous function. This Problem refers to the so-called _primal_ formulation of the minimization problem and ![x in R^N](https://latex.codecogs.com/svg.latex?x\in\mathbb{R}^N) is known as the primal variable we are interested in recovering.


[generalPrimalFormulation]: https://latex.codecogs.com/svg.latex?\min_{x}&space;F(Ax) "Problem formulation"
[matA]: https://latex.codecogs.com/svg.latex?A "A"
[funcF]: https://latex.codecogs.com/svg.latex?F" "F"

## Authors
* Hendrik Dirks ([hendrik.dirks@wwu.de](mailto:hendrik.dirks@wwu.de))*
* Lars Haalck ([lars.haalck@wwu.de](mailto:lars.haalck@wwu.de))*

\*Institute for Computational and Applied Mathematics
University of Muenster, Germany

## License
**FlexBox** is copyright ©2016-2017 by Hendrik Dirks.
If you plan to distribute the software (commercially or not), please contact Hendrik Dirks for more information.

## Dependencies
In order to use the MATLAB version of FlexBox the following requirements should be met:
* Matlab >= R2015b
* Image Processing Toolbox

In order to use the C++/CUDA version please refer to the repository: https://github.com/HendrikMuenster/flexBox_CPP

## Usage
We recommend to look at the provided examples in the folder examples/.

## C++/CUDA Module

**FlexBox** comes with a C++ module which can be used stand-alone or together with MATLAB via MEX-interfaces. The C++ module can be found at https://github.com/HendrikMuenster/flexBox_CPP and it is included in this repository as a submodule in the directory flexBox_CPP. For installation instructions please read the README in the linked repository.

## Citation

If you use this toolbox please use the following citation
```
@Article{dirks2015flexbox,
  Title         = {A Flexible Primal-Dual Toolbox},
  Author        = {Dirks, Hendrik},
  Journal       = {ArXiv e-prints},
  Year          = {2016},
  Month         = mar,
  Keywords      = {Mathematics - Optimization and Control, Computer Science - Computer Vision and Pattern Recognition, Computer Science - Mathematical Software, I.4, G.1.6, G.4},
  Primaryclass  = {math.OC}
}
```
A preprint of the article can be found at http://arxiv.org/abs/1603.05835

## Reporting Bugs
In case you experience any problems, please create an issue at https://github.com/HendrikMuenster/flexBox/issues or
https://github.com/HendrikMuenster/flexBox_CPP/issues
