## Code to acompany: [Simple and maximally robust processes with no classical common-cause or direct-cause explanation](https://arxiv.org/abs/2101.11630)

#### Marcello Nery, Marco Túlio Quintino, Philippe Allard Guérin, Thiago O. Maciel, Reinaldo O. Vianna

This repository contains all the code used to calculate the numerical results presented in the article "*Simple and maximally robust processes with no classical common-cause or direct-cause explanation*", Marcello Nery, Marco Túlio Quintino, Philippe Allard Guérin, Thiago O. Maciel, Reinaldo O. Vianna, [	arXiv:2101.11630](https://arxiv.org/abs/2101.11630)".

### MATLAB code requires:
- [CVX](http://cvxr.com/cvx/download/) - A free MATLAB toolbox for rapid prototyping of optimization problems
- [QETLAB](http://www.qetlab.com/) - A free MATLAB toolbox for quantum entanglement theory
### The CVX was tested with the solvers:
- [Mosek](http://docs.mosek.com/9.0/toolbox/index.html) - The MOSEK optimization toolbox for MATLAB manual
- [SeDuMi](http://sedumi.ie.lehigh.edu/?page_id=58) - A Matlab toolbox for optimization over symmetric cones
- [SDPT3](https://www.math.cmu.edu/~reha/sdpt3.html) - A Matlab software package for semidefinite programming

For instance, installing CVX also installs the solvers (Mosek, SeDuMi, SDPT3). 

For a fast reproduction of the results, follow the following steps:

1. Add the whole repository to the MATLAB path (including the dependencies);
2. Execute the command ```results_summary``` in the MATLAB terminal. 

This command will calculate every robustness for every process (with low dimensions) presented in the work. For calculating the robustnesses of processes with high dimensions (![equation](http://latex.codecogs.com/gif.latex?W_{\textrm{FB}}), ![equation](http://latex.codecogs.com/gif.latex?W_{339})), execute ```results_summary_high_dims``` in the MATLAB terminal. This script demands a high amount of memory from the computer, so use this comand carefully.

#### Code

- [results_summary.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/results_summary.m):
Script that calculates every robustness of every process with small dimensions presented in our work (![equation](http://latex.codecogs.com/gif.latex?W_{222}), ![equation](http://latex.codecogs.com/gif.latex?W_{224}), ![equation](http://latex.codecogs.com/gif.latex?W_{\textrm{PPT}}), ![equation](http://latex.codecogs.com/gif.latex?W_{\textrm{SEP}}), ![equation](http://latex.codecogs.com/gif.latex?W_{\textrm{MRSR}}) ![equation](http://latex.codecogs.com/gif.latex?W_{\textrm{FB}}), ![equation](http://latex.codecogs.com/gif.latex?W_{339}))


- [bipartite_ordered_test.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/bipartite_ordered_test.m): 
This function tests if a given matrix represents a valid bipartite ordered process.

- [bra.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/bra.m):
Function that generates a bra vector in the computational basis.

- [ccdc_robustness_inner.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ccdc_robustness_inner.m):
Function that returns the robustness (either generalized or white noise) and optimal non-classical CCDC witness of a given bipartite ordered process using the inner approximation.

- [ccdc_robustness_outer.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ccdc_robustness_outer.m):
Function that returns the robustness (either generalized or white noise) and optimal non-classical CCDC witness of a given bipartite ordered process using the outer approximation.

- [ccdc_robustness_summary.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ccdc_robustness_summary.m):
Function that calculates the robusnesses of a given process with different methods (primal and dual, inner and outer approximation).

- [ccdc_seesaw.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ccdc_seesaw.m):
Function that numerically searches for a bipartite ordered process that attains the maximum (generalized or white noise) robustness among the processes with the same dimensions.

- [ccdc_tripartite_robustness_inner.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ccdc_tripartite_robustness_inner.m)
Function that returns the robustness (either generalized or white noise) and optimal non-classical CCDC witness of a given tripartite ordered process using the inner approximation.

- [ccdc_tripartite_robustness_outer.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ccdc_tripartite_robustness_outer.m)
Function that returns the robustness (either generalized or white noise) and optimal non-classical CCDC witness of a given tripartite ordered process using the outer approximation.

- [ket.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ket.m):
Function that generates a ket vector in the computational basis.

- [ketbra.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/ketbra.m):
Function that generates a ketbra operator in the computational basis.

- [non_ccdc_process.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/non_ccdc_process.m):
Function that generates one of the processes presented in our work.

- [non_ccdc_processes_generator.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/non_ccdc_processes_generator.m):
Function that generates all processes presented in our work.

- [process_max_wit_viol.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/process_max_wit_viol.m):
Function that searches for a process that maximally violates a given non-classical CCDC witness (used in the see-saw).

- [random_bipartite_ordered_process.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/random_bipartite_ordered_process.m):
Function that randomly generates a valid bipartite ordered process.

- [NRandomPureStates.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/NRandomPureStates.m)
Function that returns a set of random pure states of a given dimension (used in every code that involves inner approximation).

- [traceandrep.m](https://github.com/marcellongvb/non_ccdc_processes/blob/master/code/traceandrep.m)
Function that executes the operation of partial trace over a given set of subsystems and replacing it by normalized identities.
