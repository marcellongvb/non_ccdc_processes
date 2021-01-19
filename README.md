## Code to acompany: [Simple and maximally robust processes with no classical common-cause or direct-cause explanation]{https://arxiv.org/xxx}

#### Marcello Nery, Marco Túlio Quintino, Philippe Allard Guérin, Thiago O. Maciel, Reinaldo O. Vianna

This repository contains all the codes used to calculate the numerical results presented in the article "*Simple and maximally robust processes with no classical common-cause or direct-cause explanation*", Marcello Nery, Marco Túlio Quintino, Philippe Allard Guérin, Thiago O. Maciel, Reinaldo O. Vianna, [arXiv:XXX](https://arxiv.org/xxx)".

MATLAB code requires:
- [Yalmip]{https://yalmip.github.io/} - A free toolbox for modeling and optimization in MATLAB
- [CVX]{http://cvxr.com/cvx/download/) - A free MATLAB toolbox for rapid prototyping of optimization problems
- [Mosek]{http://docs.mosek.com/9.0/toolbox/index.html} - The MOSEK optimization toolbox for MATLAB manual
- [SeDuMi]{http://sedumi.ie.lehigh.edu/?page_id=58} - A Matlab toolbox for optimization over symmetric cones
- [SDPT3]{https://www.math.cmu.edu/~reha/sdpt3.html} - A Matlab software package for semidefinite programming
- [QETLAB]{http://www.qetlab.com/} - A free MATLAB toolbox for quantum entanglement theory

For instance, installing CVX also installs the solvers (Mosek, SeDuMi, SDPT3). 

For a fast reproduction of the results, follow the following steps:

1. Add the whole repository to the MATLAB path (including the dependencies);
2. Execute the command ```results_summary``` in the MATLAB terminal. 

This command will calculate every robustness for every process (with low dimensions) presented in the work. For calculating the robustnesses of processes with high dimensions ($W_{\textrm{FB}}$, $W_{339}$), execute ```results_summary_high_dims``` in the MATLAB terminal. This script demands a high amount of memory from the computer, so use this comand carefully.

The MATLAB codes of this repository are: 


