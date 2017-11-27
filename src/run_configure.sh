#!/bin/sh

# Executable must be 32-bit.
# -std=c99  is necessary to get rid of long long warning in netcdf.h
#    But that causes problems with other stuff

# Make sure this is for the x86 (not x86_64) version of NetCdf
# Assumes nc-config is in the path.
#NC_CONFIG = "nc-config"
# Assumes nc-config is not in the path.
NC_CONFIG="/usr/local/netcdf/x86-gcc/netcdf-current/bin/nc-config"

# Get the commands that were used to compile netcdf
CC=`${NC_CONFIG} --cc`
LDFLAGS=`${NC_CONFIG} --libs`
CFLAGS=`${NC_CONFIG} --cflags`

# No NetCDF Available
#CFLAGS="${CFLAGS} -m32 -static"
# Have NetCDF and HDF
CFLAGS="${CFLAGS} -m32 -static -DHAVE_NETCDF_SUPPORT -DHAVE_HDF_SUPPORT"

# This should work, but LDFLAGS is not handled properly
#echo "CC=${CC}"
#echo "LDFLAGS=${LDFLAGS}"
#echo "CFLAGS=${CFLAGS}"
#./configure CC="${CC}" CFLAGS="${CFLAGS}" LDFLAGS="${LDFLAGS}"

echo "CC=${CC}"
echo "CFLAGS=${CFLAGS} ${LDFLAGS}"
./configure CC="${CC}" CFLAGS="${CFLAGS} ${LDFLAGS}"
