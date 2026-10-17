"""Write sag_mapping_summary.tsv from finished PAF files (classify_paf.py logic)."""
import glob,os,csv
exec(open('classify_paf.py').read().split("if __name__")[0])
runs={r['run_accession']:r for r in csv.DictReader(open('runs.tsv'),delimiter='\t')}
o=open('sag_mapping_summary.tsv','w')
o.write('run\tsample\treads_screened\treads_mycobacterial_ge95\tFSM_reads\tFSM_reads_ge99\tnonFSM_reads\tmixed\tFSM_windows\tmax_windows_any_genome\ttop_genomes\n')
T=collections.Counter(); maxw=0
for d in sorted(glob.glob('paf/*.done')):
    f=d[:-5]; rid=os.path.basename(f)[:-4]
    c,G,W,FW=run(f)
    mw=max([len(v) for v in W.values()],default=0); maxw=max(maxw,mw)
    tot=c['FSM']+c['nonFSM']+c['mixed']
    o.write(f"{rid}\t{runs[rid]['sample_alias']}\t10000000\t{tot}\t{c['FSM']}\t{c['FSM_ge99']}\t{c['nonFSM']}\t{c['mixed']}\t{len(FW['FSMclade'])}\t{mw}\t{';'.join(f'{g}:{n}r/{len(W[g])}w' for g,n in G.most_common(4))}\n")
    T.update(c); T['runs']+=1
print(dict(T),'max windows any genome',maxw)
