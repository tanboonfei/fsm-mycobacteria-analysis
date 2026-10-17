set -u
ALN=/path/to/project/analysis_3way/three_way.codon.fa
D=/path/to/project/publication/supplementary_text/relax
O=$D/new_runs
run () {  # name tree srv
  /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=8 relax \
     --alignment $ALN --tree $2 --test Test --reference Reference \
     --code Universal --models Minimal --srv $3 \
     --output $O/$1.json > $O/$1.log 2>&1
  echo "$1 exit=$? $(date -u +%H:%M:%S)" >> $O/status.txt
}
run srv_FSM_tips_stem   $D/common_reference_FSM_tips_stem.nwk Yes &
run srv_LEP_tips_stem   $D/common_reference_LEP_tips_stem.nwk Yes &
run joint_FSM_vs_LEP    $D/joint_FSM_vs_LEP.nwk               No  &
wait
echo "ALL DONE $(date -u +%H:%M:%S)" >> $O/status.txt
