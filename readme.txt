Welcome to the RexOut program
Arthur Taylor
8/20/2009

The purpose of this program is to read a SLOSH rexfile, and extract a time
series of the surge according to the SLOSH model for a particular grid cell.

To do so, it uses a .pts file which is plain text.  It contains the basin
abbreviation which is cross-checked with the rexfile to make sure that the
rexfile was run on the correct basin.  Next the .pts file contain a set of
x y pairs, separated by a space, one pair per line.

Usage: .\bin\rexout.exe [OPTION]...
Options:
-help           (1 arg) usage or help command
-V              (1 arg) version command
-rex            rex file to read from
-pnt            point file to read from
-out            output file (defaults to stdout)
-bsnDir         directory containing sloshbsn info
-basin          SLOSH basin to use with this rexfile.
-style          style of output
                0 => original output format
                1 => original output but with FORTRAN fixed fields
                2 => fixed field, wind output, 1 point per section
                     Note: 2 requires the -bsnDir option
                3 => fixed field, wind output, model depth, 1 pnt per section
                     Note: 3 requires the -bsnDir option and -basin options
-incEnv         include the envelope frame

Simplest way to run requires: -rex, -pnt, -style
