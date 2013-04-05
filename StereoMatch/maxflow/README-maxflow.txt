This directory contains Vladimir Kolmogorov's <vnk@cs.cornell.edu>
interface to Andrew Goldberg's push-relabel maximum flow code.
This interface consists of maxflow.h and maxflow.c, as well as
a modified version of parser_fl.c called parser_flg.c.

To enable this code to compile (and hence the graph cut algorithm to run)
you must first download the push-relabel flow code from:
    http://www.intertrust.com/star/goldberg/soft.html
(read the Copyright and then download the PRF (tar format) code).
Alternatively, you can simply not define PRF_DOWNLOADED in maxflow.c,
and not run the graph cut algorithm.

From prf.tar, extract the following files into the same directory as maxflow.c:

05/30/2000  06:08p                 970 types_pr.h
05/30/2000  06:08p               3,193 phase2.c
05/30/2000  06:08p               9,729 h_prf.c

You can also extract the other *_prf.c and df.c files if you wish to
play with these other algorithm variants.

After you have done this, define the PRF_DOWNLOADED macro parameter,
and recompile maxflow.c.