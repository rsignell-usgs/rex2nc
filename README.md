rex2nc
======

Justin Davis (@JRDavisUF) slosh rex file to netcdf converter

Simply using "./configure" will NOT compile the code using the netcdf routines in the code. There is nothing in ANY of the configure stuff that I see to provide the needed CPP variable "HAVE_NETCDF_SUPPORT" 

If one wants NetCDF support, one must supply the needed NetCdf configuration information to the configure script (as what NETCDF=enable would normally do) derived from the 32-bit version of nc-config. 

I put the needed info in [another script](./src/run_configure.sh) which specifies variables and then runs configure.  Not how I specifically defined the HAVE_NETCDF_SUPPORT CPP variable to ensure the correct parts of the code get compiled.

Don't forget, one must build a 32-bit executable...which means one must use the nc-config program from a 32-bit install of NetCDF.

jrd
