"""Opsin neighbourhood synteny from the uniform Bakta annotation of +-15 kb windows."""
import pandas as pd, itertools, numpy as np, collections
O='/path/to/project/analyses/extended/B_opsin'
b=pd.read_csv(O+'/synteny/bakta_windows/windows.tsv',sep='\t',comment='#',header=None,names=['seq','type','start','stop','strand','tag','gene','product','xref'])
b=b[b.type=='cds'].sort_values(['seq','start'])
cl=pd.read_csv(O+'/synteny/window_protein_clusters.tsv',sep='\t',header=None,names=['rep','tag']); c2=dict(zip(cl.tag,cl.rep))
L=pd.read_csv(O+'/screen/opsin_loci_lineages.tsv',sep='\t'); L=L[L.bakta_tag.notna()]
L['grp']=L.apply(lambda x:f"{'FSM' if x.clade=='FSM' else 'nonFSM'}-{x.lineage}",axis=1)
K=5; rows=[]; neigh={}
isblh=b['product'].str.contains('carotene 15,15',na=False)
for _,x in L.iterrows():
    g=b[b.seq==x.locus].reset_index(drop=True)
    if x.bakta_tag=='none':
        i=None
    else:
        i=g.index[g.tag==x.bakta_tag][0]
    if i is None: continue
    st=g.loc[i,'strand']
    win=g.iloc[max(0,i-K):i+K+1]
    rel={}
    for j,r in win.iterrows():
        if j==i: continue
        off=(j-i) if st=='+' else (i-j)   # + = downstream in opsin orientation
        rel[off]=r
    neigh[x.locus]=set(c2[r.tag] for r in rel.values())
    blh=[o for o,r in rel.items() if isblh[r.name] if False] if False else [o for o,r in rel.items() if 'carotene 15,15' in str(r['product'])]
    left=sum(1 for o in rel if o<0); right=sum(1 for o in rel if o>0)
    rows.append(dict(locus=x.locus,genome=x.genome,clade=x.clade,lineage=x.lineage,sublineage=x.sublineage,grp=x.grp,
        genes_upstream=left,genes_downstream=right,complete_5_each_side=(left==K and right==K),
        blh_adjacent=int(any(abs(o)==1 for o in blh)),blh_within5=int(len(blh)>0),blh_offset=';'.join(map(str,sorted(blh))),
        **{f'g{o:+d}':str(rel[o]['product']) if o in rel else '' for o in range(-K,K+1) if o!=0}))
S=pd.DataFrame(rows); S.to_csv(O+'/synteny/opsin_neighbourhoods.tsv',sep='\t',index=False)
print(S.groupby('grp').agg(loci=('locus','size'),complete=('complete_5_each_side','sum'),blh_adjacent=('blh_adjacent','sum'),blh_within5=('blh_within5','sum')))
# immediate-neighbour products per group
for gname,gg in S.groupby('grp'):
    print('\n##',gname,len(gg))
    for o in ['g-1','g+1']:
        print(' ',o, gg[o].replace('', 'contig end').value_counts().head(5).to_dict())
# Jaccard similarity of +-5 neighbour cluster sets, within and between groups
grp=dict(zip(S.locus,S.grp)); loci=[l for l in S.locus if S.set_index('locus').loc[l,'complete_5_each_side']]
J=collections.defaultdict(list)
for a,c in itertools.combinations(loci,2):
    if S.set_index('locus').loc[a,'genome']==S.set_index('locus').loc[c,'genome']: continue
    k=tuple(sorted([grp[a],grp[c]])); u=neigh[a]|neigh[c]
    J[k].append(len(neigh[a]&neigh[c])/len(u) if u else 0)
out=pd.DataFrame([dict(group1=k[0],group2=k[1],pairs=len(v),mean_jaccard=np.mean(v),median_jaccard=np.median(v),frac_pairs_sharing_any=np.mean([x>0 for x in v])) for k,v in J.items()])
out.to_csv(O+'/synteny/neighbourhood_jaccard.tsv',sep='\t',index=False); print(out.round(3).to_string())
# clusters most often found next to opsins per group
rep2prod={}
for t,r in zip(b.tag,b['product']): rep2prod.setdefault(c2[t],r)
fr=[]
for gname,gg in S.groupby('grp'):
    cnt=collections.Counter(c for l in gg.locus for c in neigh[l])
    for c,n in cnt.most_common(8): fr.append(dict(grp=gname,cluster=c,product=rep2prod[c],loci=n,frac=n/len(gg)))
F=pd.DataFrame(fr); F.to_csv(O+'/synteny/top_neighbour_families.tsv',sep='\t',index=False); print(F.round(2).to_string())
