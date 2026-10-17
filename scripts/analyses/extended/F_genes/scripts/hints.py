"""Top TCDB and Swiss-Prot hits (DIAMOND blastp --more-sensitive) for the cluster representatives
covering >=10 FSM genomes (or >=10 sister genomes) per target; TC family names from TCDB."""
from common import *
fam={}
for l in open(SP+'refdb_dl/tcdb_families.tsv',encoding='utf-8',errors='replace'):
    p=l.rstrip('\n').split('\t')
    if len(p)>=2: fam[p[0]]=p[1].strip()
def top(fn):
    t={}
    for l in open(fn):
        f=l.rstrip('\n').split('\t')
        if f[0] not in t: t[f[0]]=f
    return t
tc=top(SP+'hint_tcdb.tsv'); sp=top(SP+'hint_sprot.tsv')
R=pd.read_csv(SP+'reps_targets.tsv',sep='\t')
R=R[(R.FSM_genomes>=10)|(R.sister_genomes>=10)]
extra=pd.DataFrame([dict(point='2_transport',gene='oppA (genefam FSM-divergent family, S14)',rep='spire_mag_00841652__k141_1370522_68',FSM_genomes=43,sister_genomes=0,eggNOG_desc='genefam MAG-panel family',PFAMs='SBP_bac_5',seed='931627.MycrhDRAFT_5175'),
                   dict(point='2_transport',gene='oppA (genefam K15580 family, S14)',rep='spire_mag_02088287__k141_1054254_1',FSM_genomes=22,sister_genomes=0,eggNOG_desc='genefam MAG-panel family',PFAMs='SBP_bac_5',seed='1032480.MLP_13720'),
                   dict(point='2_transport',gene='ptsG (genefam, S14)',rep='spire_mag_01948467__k141_38987_4',FSM_genomes=91,sister_genomes=0,eggNOG_desc='genefam; KO K02755/K02756/K02757',PFAMs='PTS_EIIA_1,PTS_EIIB,PTS_EIIC',seed='1108045.GORHZ_164_00470')])
R=pd.concat([R,extra])
rows=[]
for _,r in R.iterrows():
    t=tc.get(r.rep); s=sp.get(r.rep)
    tcid=t[1].split('|')[-1] if t else ''
    f2='.'.join(tcid.split('.')[:3]) if tcid else ''
    rows.append(dict(point=r.point,gene=r.gene,rep=r.rep,FSM_genomes=r.FSM_genomes,sister_genomes=r.sister_genomes,PFAMs=r.PFAMs,
        TCDB_top=f"{tcid} ({t[2]}% id, {int(t[3])}/{t[4]} aa) {t[8][:70]}" if t else 'no hit',TC_family=fam.get(f2,''),
        SwissProt_top=f"{s[8][:90]} ({s[2]}% id, {int(s[3])}/{s[4]} aa)" if s else 'no hit'))
df=pd.DataFrame(rows)
df.to_csv(F+'TableS16_new_transporter_substrate_hints.tsv',sep='\t',index=False)
pd.set_option('display.width',400); pd.set_option('display.max_colwidth',110)
print(df[df.point.isin(['2_transport'])|df.gene.isin(['pnuC','betP/betT','thiM','znuA','amtB','feoB','glpF','thiBPQ','bfr','mrpD'])][['gene','FSM_genomes','sister_genomes','TCDB_top','TC_family','SwissProt_top']].to_string())
