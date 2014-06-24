Welcome to the RexOut (JRD modified) program
Justin R. Davis
8/1/2011

The purpose of this modified program is to read a SLOSH rexfile, and
extract a spatial field or a time series of the surge simulated using
SLOSH in either IMEDS ASCII or NetCDF format.

Usage: ./rexout [OPTION]...

Options:
-help           (1 arg) usage or help command [MODIFIED OPTION]
-V              (1 arg) version command
-rex            rex file to read from
-pnt            point file to read from
-bnt            bnt file to read from [NEW OPTION]
-out            output file (defaults to stdout)
-bsnDir         directory containing sloshbsn info
-basin          SLOSH basin (path/filename) to use with this rexfile.
-basinName      SLOSH basin name [NEW OPTION]
-style          style of output
                0 => original output format
                1 => original output but with FORTRAN fixed fields
                2 => fixed field, wind output, 1 point per section
                     Note: 2 requires the -bsnDir option
                3 => fixed field, wind output, model depth, 1 pnt per section
                     Note: 3 requires the -bsnDir option and -basin options
                4 => NetCDF output format (time series) [NEW OPTION]
                5 => NetCDF output format (spatial)     [NEW OPTION]
                6 => IMEDS  output format for water level (time	series) [NEW OPTION]

-incEnv         include the envelope frame

Simplest way to run requires: -rex, -pnt, -style

================================================================================
Example Usage of New Options:

Style 4 (NetCDF Time Series)
./rexout -rex camilbix.rex -pnt bix.pts    -out camil.bix.ts.nc         -bsnDir . -style 4 -bnt sloshdsp.bnt

Style 5 (NetCDF Spatial)
./rexout -rex camilbix.rex -basinName hbix -out camil.bix.spatial.nc    -bsnDir . -style 5 -bnt sloshdsp.bnt

Style 6 (IMEDS ASCII Time Series)
./rexout -rex camilbix.rex -pnt bix.pts    -out camil.bix.imeds.txt     -bsnDir . -style 6 -bnt sloshdsp.bnt


bix.pts is the list of stations at which to output data.

A list of points is not required for Style 5 as all points are output...however,
  the basinName must be supplied for this style (ie, first line of .pts file).

