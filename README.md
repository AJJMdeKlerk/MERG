# MERG

MERG is a library for calculating the out of equilibrium dynamics following a
quantum quench in the interaction strength of the repulsive Lieb-Liniger model
using numerical renormalization group methods. In particular, it contains an
implementation of the High Overlap State Truncation Scheme and the Matrix
Element Renormalization Group algorithms introduced in
[https://arxiv.org/pdf/1911.11101.pdf](https://arxiv.org/pdf/1911.11101.pdf). 

If you have any questions, please do not hesitate to send me an email at
[a.j.j.m.deklerk@uva.nl](mailto:a.j.j.m.deklerk@uva.nl).

# Installing the code

1. Clone this repository by running `` $ git clone https://github.com/AJJMdeKlerk/MERG.git ``
2. Install TBB (v >= 2021.2.0), Eigen (v >=3.3.9), as well as g++-11 and make
   sure that the paths in the INCLUDES and LIBS variables in the Makefile refer
   to where you installed these external libraries. 
3. In the cloned repository, run `` $ make ``, which creates the executables in
   the build directory.

# Running the code

The executables are generated from the .cpp files in the mains folder. An
example main file is there to serve as a template. An executable can be run by
executing `` $ build/insert_name_executable `` in the base directory of the
project. 

# Word of caution 

The current version of the code does not impose restrictions on the RAM usage. 
Therefore it is advisable to keep an eye on the RAM usage when considering large
numbers of states and/or large step sizes in the NRG-routines. 

# Acknowledgements

This code was written as part of my PhD under the supervision of [J-S
Caux](https://jscaux.org/). The implementation of the Lieb-Liniger model as well
as the scanning routines are based on ideas originating from
[ABACUS](https://aip.scitation.org/doi/10.1063/1.3216474). For the development
and implementation of the numerical renormalization group routines, I
collaborated with [N.J. Robinson](https://neiljrobinson.com/), which resulted in
a [paper](https://arxiv.org/pdf/1911.11101.pdf) currently under review at
[SciPost Physics](https://scipost.org/SciPostPhys). 

Credit should also go to [Michael Crawford](https://github.com/mbcrawfo) for
providing a useful template of a [generic
makefile](https://github.com/mbcrawfo/GenericMakefile) for C++ projects on which
the Makefile in this project is based.

# License

The code in this repository is licensed under the [GNU General Public
License](https://www.gnu.org/licenses/gpl-3.0.html), version 3.
