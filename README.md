INTRODUCTION
------------

The package contains an implementation of the Covariance Matrix
Adaption Evolution Strategy (CMAES). Implemented is the basic strategy with
box constraints using C++ and the Eigen package for matrix calculations.
We follow [1] for the pure CMAES extended by a "boundary handling" as decribed
in [2]. Both papers are linked on 'The CMA Evolution Strategy' Website [3]


USAGE
-----

Generating Executables
----------------------

The package consists of the files
- README.md             // this file
- cmaes.cpp             // the main source of the CMAES algorithm
- auxiliaries.cpp/hpp   // some auxiliary functions required by CMAES
- eigenreq.hpp          // path to required EIGEN algebra package and some defs
- testfunctions.cpp     // collection of testfunctions for optimization
- Makefile              // generates executables
- nIter0.dat            // sample interface file (with operational settings)
- rosenbrock.dat        // example instance file to be used as cmaes argument

In order to generate executables, the "Eigen" algebra package must be available
on the system. It can be downloaded from
  http://eigen.tuxfamily.org
as ".tar.gz" file and simply extracted to some place on your system without any
installation process (cf. http://eigen.tuxfamily.org/dox/GettingStarted.html).
The Makefile must be adapted to contain the right path then, i.e., you have to
change the Makefile line spelling
  EIGEN = -I /gfs1/work/shmvolki/eigen-eigen-36fd1ba04c12/Eigen/
to contain your own path to the Eigen folder.
The "cmaes" executable should then be generated by typing "make".

Problem Instances
-----------------

Problem instances are provided as a command line argument to cmaes.
An example is given by "rosenbrock.dat". The file is structured as follows
line 1: number n of function parameters to be optimized
line 2: population size (number of sampled parameter sets in each iteration)
line 3: maximum number of iterations
line 4: seed for the random number generator used by cmaes
line 5: name of the objective function (implemented in testfunctions.cpp)
line 6 to 6+n-1: lower and upper bound for each parameter (two numbers per row)
line 6+n to 6+2n-1: optional indicators (0 or 1) for stochastic parameters


LITERATURE
----------

[1] N. Hansen (2016). The CMA Evolution Strategy: A Tutorial. arXiv.

[2] N. Hansen, A.S.P. Niederberger, L. Guzzella and P. Koumoutsakos (2009).
    A Method for Handling Uncertainty in Evolutionary Optimization with an
    Application to Feedback Control of Combustion.
    IEEE Transactions on Evolutionary Computation, 13(1), pp. 180-197

[3] CMAES website: https://www.lri.fr/~hansen/cmaesintro.html
