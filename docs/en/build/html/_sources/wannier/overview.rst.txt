Overview
========

In this document, we introduce how we compute downfolded models
with mVMC or HPhi++ in conjunction to
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_.

.. math::

   \begin{aligned}
   {\cal H} &=
   \sum_{R, R', i, j, \sigma}
   \left(t_{(R'-R) i j} - t_{(R'-R) i j}^{\rm DC}\right)
   c_{R' j \sigma}^{\dagger} c_{R i \sigma}
   \nonumber \\
   &+ \sum_{R, i}
   U_{0 i j} n_{R i \uparrow} n_{R i \downarrow}
   + \sum_{(R, i) < (R', j)}
   U_{(R'-R) i j} n_{R i} n_{R' j}
   - \sum_{(R, i) < (R', j)}
   J_{(R'-R) i j} (n_{R i \uparrow} n_{R' j \uparrow}
   + n_{R i \downarrow} n_{R' j \downarrow})
   \nonumber \\
   &+ \sum_{(R, i) < (R', j)}
   J_{(R'-R) i j} (
   c_{R i \uparrow}^{\dagger} c_{R' j \downarrow}^{\dagger}
   c_{R i \downarrow} c_{R' j \uparrow} +
   c_{R' j \uparrow}^{\dagger} c_{R i \downarrow}^{\dagger}
   c_{R' j \downarrow} c_{R i \uparrow} )
   \nonumber \\
   &+ \sum_{(R, i) < (R', j)}
   J_{(R'-R) i j} (
   c_{R i \uparrow}^{\dagger} c_{R i \downarrow}^{\dagger}
   c_{R' j \downarrow} c_{R' j \uparrow} +
   c_{R' j \uparrow}^{\dagger} c_{R' j \downarrow}^{\dagger}
   c_{R i \downarrow} c_{R i \uparrow} ),
   \\
   t_{0 i i}^{\rm DC} &\equiv \frac{1}{2}U_{0 i i} D_{0 i i}
   + \sum_{(R, j) (\neq 0, i)} U_{R i j} D_{0 j j}
   - \frac{1}{2} \sum_{(R, j) (\neq 0, i)} J_{R i j} D_{0 j j},
   \\
   t_{R i j}^{\rm DC} &\equiv \frac{1}{2} J_{R i j} (D_{R i j} + 2 {\rm Re} [D_{R i j}])
   -\frac{1}{2}  U_{R i j} D_{R i j},
   \quad (R, j) \neq (0, i),
   \\
   D_{R i j} &\equiv \sum_{\sigma}
   \left\langle c_{R j \sigma}^{\dagger} c_{0 i \sigma}\right\rangle_{\rm KS}.
   \end{aligned}

Prerequisite
------------

We compute the Kohn-Sham orbitals with
`QuantumESPRESSO <http://www.quantum-espresso.org/>`_
or
`xTAPP <http://xtapp.cp.is.s.u-tokyo.ac.jp/>`_,
and obtain the Wannier function, the dielectric function,
the effective interaction with RESPACK,
and simulate quantum lattice models with
mVMC or HPhi++.
Therefore, these programs must be available in our machine.
