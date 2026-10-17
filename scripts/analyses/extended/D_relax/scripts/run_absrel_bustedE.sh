#!/bin/bash
# aBSREL on FSM stem (IntNode024) and leprosy stem (IntNode017); BUSTED-E (error sink) on all branches
D=/path/to/project/analyses/extended/D_relax
cd $D/absrel && OMP_NUM_THREADS=2 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=2 absrel --alignment $D/inputs/three_way.codon.fa --tree $D/inputs/stems_test.nwk --branches Test --code Universal --srv No --multiple-hits None --output $D/absrel/absrel_stems.json > $D/absrel/absrel_stems.log 2>&1 &
cd $D/bustedE && OMP_NUM_THREADS=3 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=3 busted --alignment $D/inputs/three_way.codon.fa --tree $D/inputs/plain.nwk --branches All --code Universal --srv No --error-sink Yes --output $D/bustedE/bustedE_all.json > $D/bustedE/bustedE_all.log 2>&1 &
wait
