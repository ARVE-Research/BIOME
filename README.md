# BIOME
This is a catch-all repository for ongoing development in the BIOME series of models.

The following software needs to be available in the path when building:

- Autotools (including Autoconf, Automake, M4, and libtool)
- pkg-config
- a Fortran compiler
- HDF5 (only C library required)
- netCDF compiled against the above HDF5, including C and Fortran libraries

On first install, set up configure with `autoreconf -if`. Then `./configure` and `make`.
