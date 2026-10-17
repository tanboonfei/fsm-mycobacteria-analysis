#!/bin/bash
# Prodigal (2.6.3, -p single, as run_prodigal.sh) for genomes without a stored call,
# then miniprot (0.18-r281) of the core reference proteins against every genome.
set -u
D=/path/to/project/analyses/extended/C_pseudogene
MP=/path/to/conda/envs/miniprot/bin/miniprot
cd $D
# genome<TAB>fasta<TAB>prodigal
tail -n +2 results/manifest.tsv | awk -F'\t' '{print $1"\t"$6"\t"$7}' > work/jobs.tsv
awk -F'\t' '$3==""{print $1"\t"$2}' work/jobs.tsv | xargs -P 6 -L 1 bash -c '
  out='$D'/work/prodigal/$0.gff.gz; [ -s "$out" ] && exit 0
  case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac | prodigal -q -p single -f gff -o /dev/stdout 2>/dev/null | gzip > "$out.tmp" && mv "$out.tmp" "$out"'
while IFS=$'\t' read -r g fa pr; do
  out=work/mp/$g.paf.gz; [ -s "$out" ] && continue
  $MP -t 6 --gff -N 3 "$fa" work/core_ref_proteins.faa 2>>work/mp/miniprot.log \
    | grep -E "^##PAF|	mRNA	" | sed 's/\tcs:Z:[^\t]*//' | gzip > $out.tmp && mv $out.tmp $out
done < work/jobs.tsv
echo done
