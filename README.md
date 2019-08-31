
# nuTRLan version 0.2

The up-to-date source code of this package is available at
* <https://codeforge.lbl.gov/projects/trlan/>.

Additional information about the package can be found on the web at
* <http://lbl.gov/~kwu/trlan/>  and
* <http://www.nersc.gov/research/SIMON/trlan.html>.

In particular, the user's guide is available at
* <https://codeforge.lbl.gov/docman/view.php/43/24/nutrlan-ug.pdf>
* <http://lbl.gov/~kwu/ps/trlan_ug.html>.

Further problems/questions, contact
* Ichitaro Yamazaki (IYamazaki@lbl.gov),
* Kesheng John Wu (John.Wu@ACM.org), or
* Horst D Simon (hdsimon@lbl.gob).

## Orthogonalization

Lanczos vectors may not necessarily remain in the Krylov space
due to floating-point roundoff.  In cases where this is a problem,
the user can supply an additional re-orthogonalization routine.

## Installation guide

You can use the Makefile in this directory to generate nuTRLan library.
The make file flags are in file Make.inc.  Check the file to make sure
that the options are appropriate for your particular machine.

 - You will need a C compiler to compile the main body of the library.
 - You will also need interface to LAPACK and BLAS.

There are some examples of Make.inc files

*    `Make.csia` (HP Itanium machine)
*    `Make.dmx` (Linux)
*    `Make.seaborg` (Nersc IBM SP)

`make lib` will build a sequential version of the library,

`make plib` will build a MPI version of the library.

nuTRLan calls some BLAS routines. These routines are contained in CBLAS
directory.  If BLAS is not installed on the target machine, the user can
use, these unoptimized routines by  `make pblib`.

The `examples` directory contains a small set of examples.

* Example on how to compute eigenpairs of a symmetric matrix on a
   serial machine

    `simple.c` and `test_ctrlan.c`

* Example on how to compute eigenpairs of a symmetric matrix on a
   distributed memory machine

    `psimple.c` and `test_ptrlan.c`

* Example on how to compute eigenpairs of a Hermitian matrix on a
   serial machine

    `test_ztrlan.c`

* Example on how to compute eigenpairs of a Hermitian matrix on a
   distributed memory machine

    `test_pztrlan.c`

* Example on a Fortran interface to compute eigenpairs of a symmetric
   matrix on a serial machine

    `simple.f`

* Example on a Fortran interface to compute eigenpairs of a symmetric
   matrix on a distributed memory machine

    `psimple.f`

* Example on a Fortran interface to compute eigenpairs of a Hermitian
   matrix on a serial machine

    `zsimple.f`

* Example on a Fortran interface to to compute eigenpairs of Hermitian
   matrix on a parallel machine

    `pzsimple.f`
