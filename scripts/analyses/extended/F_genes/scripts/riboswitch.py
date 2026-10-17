"""Associate Rfam TPP riboswitch hits (cmsearch --cut_ga, RF00059) with the downstream tblastn locus
(same contig and strand, riboswitch 3' end within 400 bp upstream of the locus start)."""
from common import *
import collections
loci=pickle.load(open(SP+'tbn/loci.pkl','rb'))
s1=load_meta()
def rs(fn):
    R=[]
    for l in open(fn):
        if l[0]=='#': continue
        f=l.split(); a,b=int(f[7]),int(f[8]); R.append((f[0],f[9],a,b,float(f[14])))
    return R
rows=[]; per=collections.defaultdict(lambda:collections.defaultdict(set)); nrs={}
for tag,fn in (('FSM',SP+'tbn/tpp_fsm.tbl'),('RGM_sister',SP+'tbn/tpp_sis.tbl')):
    R=rs(fn); nrs[tag]=(len(R),len(set(r[0].split('|')[0] for r in R)))
    L=collections.defaultdict(list)
    for g,qs in loci[tag].items():
        for q,h in qs.items():
            if h['pid']>=30 and h['qcov']>=40: L[(h['contig'],h['strand'])].append((q,h))
    for c,st,a,b,sc in R:
        g=c.split('|')[0]; hit='none'
        for q,h in L.get((c,st),[]):
            if st=='+' and 0<=h['lo']-b<=400: hit=q
            if st=='-' and 0<=b-h['hi']<=400: hit=q   # minus strand: cm end b < start a
        rows.append(dict(set=tag,genome=g,contig=c,strand=st,start=a,end=b,score=sc,downstream_query=hit))
        per[tag][hit].add(g)
df=pd.DataFrame(rows); df.to_csv(F+'TPP_riboswitch_hits.tsv',sep='\t',index=False)
out=[]
for tag in per:
    tot=len(members(s1,tag))
    for q,gs in sorted(per[tag].items(), key=lambda x:-len(x[1])):
        withgene=sum(1 for g in members(s1,tag) if q in loci[tag].get(g,{}) and loci[tag][g][q]['pid']>=40 and loci[tag][g][q]['qcov']>=50) if q!='none' else ''
        out.append(dict(set=tag,downstream_gene=q,genomes_with_TPP_upstream=len(gs),genomes_with_gene_credible=withgene,genomes_total=tot))
pd.DataFrame(out).to_csv(F+'TPP_riboswitch_summary.tsv',sep='\t',index=False)
print(nrs); print(pd.DataFrame(out).to_string())
