"""Genome-wide Blh (tBLASTn; 28 mycobacterial Blh from opsin windows + Q9HNE6 + Q4PNI0; e<=1e-10, >=50% query span) and its position relative to opsin loci."""
import pandas as pd, glob, os
O='/path/to/project/analyses/extended/B_opsin'
s1=pd.read_csv('/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv',sep='\t'); s1['key']=s1.tip.str.replace(r'\.fa$','',regex=True)
cl=dict(zip(s1.key,s1.clade.str.split().str[0])); org=dict(zip(s1.key,s1.origin))
rows=[]
for f in glob.glob(O+'/synteny/blh_hits/*.tsv'):
    if os.path.getsize(f)==0: continue
    d=pd.read_csv(f,sep='\t',header=None,names='q s pid len qs qe ss se ev bits slen qlen'.split()); d['g']=os.path.basename(f)[:-4]; rows.append(d)
H=pd.concat(rows); H=H[(H.qe-H.qs+1)>=0.5*H.qlen]
H['lo']=H[['ss','se']].min(axis=1); H['hi']=H[['ss','se']].max(axis=1)
best=H.sort_values('bits',ascending=False).drop_duplicates(['g','s','lo'])
# collapse to loci per contig
loc=[]
for (g,c),x in best.groupby(['g','s']):
    x=x.sort_values('lo'); cur=None
    for _,r in x.iterrows():
        if cur is not None and r.lo<=cur['hi']+200: cur['hi']=max(cur['hi'],r.hi); cur['bits']=max(cur['bits'],r.bits); continue
        if cur is not None: loc.append(cur)
        cur=dict(genome=g,contig=c,lo=r.lo,hi=r.hi,bits=r.bits,best_query=r.q,pid=r.pid)
    loc.append(cur)
B=pd.DataFrame(loc)
OL=pd.read_csv(O+'/screen/opsin_loci_lineages.tsv',sep='\t'); OL=OL[OL.contig.notna()]
def rel(r):
    o=OL[(OL.genome==r.genome)&(OL.contig==r.contig)]
    if not len(o): return pd.Series(dict(on_opsin_contig=False,dist_to_opsin_bp=None,opsin_lineage=None))
    dist=[max(0,max(r.lo,a)-min(r.hi,b)) for a,b in zip(o.locus_start,o.locus_end)]
    i=min(range(len(dist)),key=lambda k:dist[k])
    return pd.Series(dict(on_opsin_contig=True,dist_to_opsin_bp=dist[i],opsin_lineage=o.lineage.iloc[i]))
B=pd.concat([B,B.apply(rel,axis=1)],axis=1)
B['clade']=B.genome.map(cl); B['origin']=B.genome.map(org)
B.to_csv(O+'/synteny/blh_genomewide_loci.tsv',sep='\t',index=False)
G=pd.DataFrame({'genome':s1.key,'clade':s1.key.map(cl),'origin':s1.key.map(org)})
G['blh']=G.genome.isin(B.genome); G['blh_near_opsin']=G.genome.isin(B[B.on_opsin_contig&(B.dist_to_opsin_bp<=2000)].genome)
gl=pd.read_csv(O+'/species/genome_lineages.tsv',sep='\t').set_index('genome')
for k in ['any','A','B']: G[k]=G.genome.map(gl[k]).fillna(0).astype(int)
G=G[G.clade.isin(['FSM','RGM','SGM'])]
print('Blh loci',len(B),'in',B.genome.nunique(),'genomes; best-hit identity range',B.pid.min(),B.pid.max())
print(pd.crosstab(G.clade,[G.blh],margins=True))
print('Blh within 2 kb of an opsin, by clade:'); print(G.groupby('clade').blh_near_opsin.sum())
print('FSM: opsin x Blh'); print(pd.crosstab(G[G.clade=='FSM']['any'],G[G.clade=='FSM'].blh))
print('non-FSM opsin carriers with Blh anywhere:',G[(G.clade!='FSM')&(G['any']==1)].blh.sum(),'of',(G[(G.clade!='FSM')&(G['any']==1)]).shape[0])
print('non-FSM opsin-negative with Blh:',G[(G.clade!='FSM')&(G['any']==0)].blh.sum(),'of',(G[(G.clade!='FSM')&(G['any']==0)]).shape[0])
x=B[B.clade=='FSM']; print('FSM Blh loci on opsin contig:',x.on_opsin_contig.sum(),'; FSM Blh best query:',x.best_query.value_counts().head(3).to_dict(),'identity median',x.pid.median())
G.to_csv(O+'/synteny/blh_genome_table.tsv',sep='\t',index=False)
from scipy.stats import fisher_exact
for c in ['FSM','RGM','SGM']:
    t=pd.crosstab(G[G.clade==c]['any'],G[G.clade==c].blh).reindex(index=[0,1],columns=[False,True],fill_value=0)
    print(c,'opsin x Blh Fisher',t.values.tolist(),fisher_exact(t.values))
