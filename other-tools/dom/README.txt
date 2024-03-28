Compiling:

cc -O dom.c -o dom util/aa/cones.o util/aa/stutest.o util/aa/bestrot.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o util/aa/matrix.o util/aa/siva.o util/aa/ql.o -lm

#include <alloca.h>
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/pdbprot.h"
#include "util/aa/incl/matplot.h"
#include "util/aa/incl/matrix.h"

Most of the call to the util routines seem to be for calculating inertial axes.
I can't see where so they can probably be dropped.

Running:

./dom pdb-file [spread nruns subdoms]

pdb-file can be full of CA (just CAs used)

spread = granularity level (default 15)

nruns = flag to set mode as..
        nruns          -3  -2  -1   0   1   2
        structures      N   S   N   S  NS  NS
        beta-bias       -   -   B   B   -   B
        where:
                structures: N = native, S = smoothed, NS = find the value of spread where N and S agree
                beta-bias:  apply the same algorithm to a network of sheet (pseudo) H-bonds

subdoms = continue parsing large domains

Output:

dom0.out = pdb coloured (temp factor field) by domain number
domN.out (N = 1,2,3...) = each domain (with the cut ends linked by a fake loop

Example:

1aoz.cas = three tightly packed beta domains

./dom 1aoz.cas 15 0
default spread (15) on smooth chain (0) with sheet bias finds 2 domain (same for native (1))

./dom 1aoz.cas 11 0
reduced spread (11) on smooth chain (0) with sheet bias finds 3 domains (same for native (1))
