#!/bin/bash
D=/path/to/project/analyses/extended/D_relax; C=$1; T=$2
OMP_NUM_THREADS=$T /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$T relax --alignment $D/controls/$C.codon.fa --tree $D/controls/$C.lbl.nwk \
   --test Test --reference Reference --code Universal --models Minimal --srv No --output $D/controls/$C.relax.json > $D/controls/$C.relax.log 2>&1
