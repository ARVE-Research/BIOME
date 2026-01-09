# BIOME
This is a catch-all repository for ongoing development in the BIOME series of models.

The following software needs to be available in the path when building:

- Autotools (including Autoconf, Automake, M4, and libtool)
- pkg-config
- a Fortran compiler
- HDF5 (only C library required)
- netCDF compiled against the above HDF5, including C and Fortran libraries

On first install, set up configure with `autoreconf -if`. Then `./configure` and `make`.

Copy the sample run parameters namelist to a local file and edit the pathnames to where your data files are stored:

`cp NorthAm.namelist NorthAm.local.namelist`

Sample run command:

`src/biome1 joboptions/NorthAm.local.namelist -2117500/-617500/-1392500/-112500 output/test.nc`

If only a single pixel is run – this is done by specifing only an `xcoord/ycoord` on the command line – the model will additionally output a text file containing selected daily meteorological variables with the filename `dmet_xcoord_ycoord.txt`. A sample plotting script for visualizing this timeseries is providied in `output/plotweather2.sh`.
