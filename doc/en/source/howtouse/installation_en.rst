Installation
============

You can download HPhi++ at the following location.

https://github.com/QLMS/HPhi/releases

You can obtain the HPhi++ directory by typing

``$ tar xzvf HPhi-xxx.tar.gz``

HPhi++ can be installed by using cmake.

.. tip::

 | Before using cmake for sekirei, you must type
 | ``source /home/issp/materiapps/tool/env.sh``

We can compile HPhi++ as::

 cd $HOME/build/hphi
 cmake -DCONFIG=gcc $PathTohphi
 make

To use ScaLAPACK library for full diagonalization, the cmake option ``-DUSE_SCALAPACK=ON`` is needed.
Here, we set a path to HPhi++ as ``$PathTohphi``
and to a build directory as ``$HOME/build/hphi``.
After compilation, ``src`` folder is constructed below a ``$HOME/build/hphi``
folder and we obtain an executable ``HPhi`` in ``src/`` directory.
When no MPI library exists in the system, an executable ``HPhi``
is automatically compiled without an MPI library.

In the above example,
we compile HPhi++ by using a gcc compiler.
We can select a compiler by using the following options:

 * ``sekirei`` : ISSP system-B \"sekirei\"
 * ``sekirei_acc``: ISSP system-B \"sekirei\" (for using MAGMA library)
 * ``fujitsu`` : Fujitsu compiler
 * ``intel`` : Intel compiler + Linux PC
 * ``gcc`` : GCC compiler + Linux PC.

An example of compiling HPhi++ by using the Intel compiler is shown as follows::

 mkdir ./build
 cd ./build
 cmake -DCONFIG=intel ../
 make

After compilation,
``src`` folder is created below the ``build`` folder and
an execute HPhi++ in the ``src`` folder.
Please note that we must delete the ``build`` folder and
repeat the above operations when we change the compiler.
