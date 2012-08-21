Fortran Utilities
=================

Various Fortran utilities, that can be included into any Fortran
program.

The modules are mostly independent of each other. Simply copy any modules that
you need into your project. Tests are in the ``tests`` directory, you can look
there for examples of usage.

License
-------

All code is MIT licensed.

Functionality
-------------

Main features:

* Types (``dp``)
* Constants (``pi``, ``e_``, ``i_``)
* Sorting
* Saving/loading 2D arrays (``savetxt``, ``loadtxt``)
* Meshes (exponential, uniform)
* Cubic splines
* Saving/loading PPM images
* Lapack interface (and a few simple f90 wrappers like ``eigh``, ``inv``)
* HDF5 interface

Requirements
------------

The modules utils and ppm in ``utils.f90`` and ``ppm.f90`` use the
``newunit`` option to ``open()``. This option is part of Fortran 2008 and
requires at least gfortran 4.5 to compile.

Contributors
------------

See the `AUTHORS
<https://github.com/certik/fortran-utils/blob/master/AUTHORS>`_ file.
