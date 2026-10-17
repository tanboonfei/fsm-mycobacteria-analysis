#!/bin/bash
# usage: single_omega.sh <FSMstem|LEPstem> <tree with {Test} on one branch> <threads>
# SingleOmega.bf (HyPhy 2.5.59) has no keyword arguments for code/alignment/tree; answers are given positionally on stdin
D=/path/to/project/analyses/extended/D_relax
printf "1\n$D/ci/$1.codon.fa\n$2\n5\n" | OMP_NUM_THREADS=$3 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$3 /usr/share/hyphy/TemplateBatchFiles/SelectionAnalyses/SingleOmega.bf > $D/ci/singleomega_$1.log 2>&1
