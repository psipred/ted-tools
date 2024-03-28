This is the "dom" program written by Dr Willie Taylor (willietaylorwork@gmail.com) and described in the following article:
William R. Taylor, Protein structural domain identification, Protein Engineering, Design and Selection, Volume 12, Issue 3, March 1999, Pages 203–216, https://doi.org/10.1093/protein/12.3.203


To compile just use make.


Running:

./dom pdb-file [spread nruns subdoms]

pdb-file can contain just C-alpha coordinates

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
