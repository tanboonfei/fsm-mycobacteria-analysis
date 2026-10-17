import re, csv, statistics as st, collections
from groups import load
G=load()
KO={l.split('\t')[0]:set(l.rstrip('\n').split('\t')[1].split(',')) for l in open('genome_ko.tsv')}
DEF={l.split('\t')[0]:l.rstrip('\n').split('\t')[1] for l in open('moddef.tsv')}
def split_top(s,sep):
    out=[];d=0;cur=''
    for ch in s:
        if ch=='(':d+=1
        elif ch==')':d-=1
        if ch in sep and d==0:
            out.append((cur,ch)); cur=''
        else: cur+=ch
    out.append((cur,None)); return out
def val(e,ks):
    e=e.strip()
    if not e: return 1.0
    parts=[p for p,_ in split_top(e,' ') if p.strip()]
    if len(parts)>1: return sum(val(p,ks) for p in parts)/len(parts)
    alts=[p for p,_ in split_top(e,',')]
    if len(alts)>1: return max(val(a,ks) for a in alts)
    # complex: + required, - optional
    toks=re.split(r'([+-])',e) if not e.startswith('(') or split_top(e,'+-')[0][1] else None
    comps=split_top(e,'+-')
    if len(comps)>1:
        req=[];prev='+'
        for c,sep in comps:
            if c.strip() and prev=='+': req.append(val(c,ks))
            prev=sep or '+'
        return sum(req)/len(req) if req else 1.0
    if e.startswith('(') and e.endswith(')'): return val(e[1:-1],ks)
    if e=='--': return None
    return 1.0 if e in ks else 0.0
def module(m,ks):
    steps=[p for p,_ in split_top(DEF[m],' ') if p.strip() and p.strip()!='--']
    v=[val(s,ks) for s in steps]; v=[x for x in v if x is not None]
    full=sum(1 for x in v if x>=0.999)
    return full,len(v)
PROD={'Lys':['M00016','M00525','M00526','M00527'],'Thr':['M00018'],'Met':['M00017'],'Cys':['M00021','M00338','M00609'],
'Ser':['M00020'],'Val/Ile':['M00019'],'Ile(2-oxobutanoate)':['M00570','M00535'],'Leu':['M00432'],'Pro':['M00015'],
'Orn':['M00028'],'Arg(Orn=>Arg)':['M00844','M00845'],'His':['M00026'],'Chorismate':['M00022'],'Trp':['M00023'],'Phe':['M00024','M00910'],
'Tyr':['M00025','M00040'],'NAD':['M00115','M00912'],'Pantothenate':['M00119'],'CoA':['M00120'],'Biotin':['M00123','M00573','M00577','M00950'],
'Thiamine':['M00127','M00895'],'Riboflavin':['M00125'],'PLP':['M00124','M00916'],'THF':['M00126','M00840','M00841'],'Heme':['M00121','M00926'],
'Cobalamin(corrin ring)':['M00925','M00924'],'Cobalamin(final)':['M00122'],'Menaquinone':['M00116'],'Lipoate':['M00881'],'MoCo':['M00880'],
'Purine(IMP)':['M00048'],'Pyrimidine(UMP)':['M00051'],'Ergothioneine':['M01062'],'Siroheme':['M00846']}
SINGLE={'Gly(glyA)':['K00600'],'Gln(glnA)':['K01915'],'Asn(asnB/asnA)':['K01953','K01914']}
rows=[]
for g,meta in G.items():
    ks=KO.get(g,set()); r={'genome':g,'group':meta['group'],'size':meta['size'],'comp':meta['comp']}
    for p,ms in PROD.items():
        best=None
        for m in ms:
            f,n=module(m,ks); miss=n-f
            if best is None or miss<best[0]: best=(miss,n)
        r[p]=best[0]  # missing steps
    for p,k in SINGLE.items(): r[p]=0 if any(x in ks for x in k) else 1
    rows.append(r)
cols=['genome','group','size','comp']+list(PROD)+list(SINGLE)
with open('pathway_missing_steps.tsv','w') as o:
    o.write('\t'.join(cols)+'\n')
    for r in rows: o.write('\t'.join(str(r[c]) for c in cols)+'\n')
# summary: % genomes with pathway complete (0 missing) or <=1 missing, completeness>=90
order=['FSM','RGM_sister','RGM_other','SGM_env','NTM_opportunistic','Leprosy','SGM_obligate']
prods=list(PROD)+list(SINGLE)
with open('pathway_summary.tsv','w') as o:
    o.write('pathway\t'+'\t'.join(f'{g}_complete%\t{g}_le1missing%' for g in order)+'\n')
    for p in prods:
        out=[p]
        for g in order:
            sub=[r for r in rows if r['group']==g and r['comp']>=90]
            n=len(sub)
            out+= [f"{100*sum(r[p]==0 for r in sub)/n:.0f}", f"{100*sum(r[p]<=1 for r in sub)/n:.0f}"] if n else ['NA','NA']
        o.write('\t'.join(out)+'\n')
print('n per group (comp>=90):',{g:sum(1 for r in rows if r['group']==g and r['comp']>=90) for g in order})
