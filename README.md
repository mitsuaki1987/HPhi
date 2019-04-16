HPhi++
======

A numerical solver package for a wide range of quantum lattice models including Hubbard-type itinerant electron hamiltonians, quantum spin models, and Kondo-type hamiltonians for itinerant electrons coupled with quantum spins. The Lanczos algorithm for finding ground states and newly developed Lanczos-based algorithm for finite-temperature properties of these models are implemented for parallel computing (hybrid parallelization with OpenMP and MPI). A broad spectrum of users including experimental researchers is cordially welcome.

### Methods
Lanczos algorithm, thermal pure quantum state, full diagonalization  

### Target models
Hubbard model, Heisenberg model, Kondo lattice model, Kitaev model, Kitaev-Heisenberg model, multi-orbital Hubbard model

### Available physical quantities
specific heat, susceptibility, ground state energy, free energy, structure factors


## Requirement
C compiler (intel, Fujitsu, GNU, etc. )  
LAPACK library (intel MKL, Fujitsu, ATLAS, etc.)  
MPI library (If you do not use MPI, it is not necessary.)

## Install

You can install HΦ and also get a manual for HΦ from a [release note](https://github.com/issp-center-dev/HPhi/releases).

## Pre-Installed system
- System B of The Institute Solid State Physics jpint-use supercomputer (cleck [here](http://www.issp.u-tokyo.ac.jp/supercom/visitor/x92nxz/hphi) for details).

## Licence

The distribution of the program package and the source codes for HPhi follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 

## Official page
[HPhi portal site](http://issp-center-dev.github.io/HPhi/index_en.html)

## Author
Youhei Yamaji, Takahiro Misawa, Synge Todo, Kota Ido, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Naoki Kawashima.
