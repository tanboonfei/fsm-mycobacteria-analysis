"""Pick tblastn query proteins: RGM isolate (M. smegmatis GCF_001457595.1 first, then other closed RGM,
then sister RGM, then M. tuberculosis H37Rv) carrying the KO/eggNOG name; FSM-specific types from a
high-completeness FSM MAG."""
from common import *
d=pickle.load(open(SP+'sets.pkl','rb'))
s1=load_meta()
pref=['GCF_001457595.1','GCA_900078775.1','GCF_000317305.3','GCF_900637195.1','GCA_000230935.2','GCF_000195955.2']
fsmpref=['spire_mag_02088284','spire_mag_02088330','spire_mag_02762480','spire_mag_01882647','spire_mag_00139887','spire_mag_02683911']
Q=[ # (query_name, key, source)
('proP','K:K03762','rgm'),('betP_RGM','K:K05020','rgm'),('betP_FSM','K:K02168','fsm'),('treS','K:K05343','rgm'),('treY','K:K06044','rgm'),('treZ','K:K01236','rgm'),
('otsA','K:K00697','rgm'),('mrpA','K:K05565','rgm'),('kup','K:K03549','rgm'),('mscL','K:K03282','rgm'),('proX','K:K02002','rgm'),('ectC','K:K06720','rgm'),
('ptsG_FSM','N:ptsG','fsm'),('nagE','K:K02804','rgm'),('pnuC_FSM','K:K03811','fsm'),('oppA_FSM','K:K15580','fsm'),('potD','K:K11069','rgm'),
('amtB','K:K03320','rgm'),('nirB','K:K00362','rgm'),('nirD','K:K00363','rgm'),('narK','K:K02575','rgm'),('cynS','K:K01725','rgm'),('urtA','K:K11959','rgm'),('ureC','K:K01428','rgm'),('glnA','N:glnA','rgm'),
('feoB','K:K04759','rgm'),('efeU','K:K07243','rgm'),('efeO','K:K07224','rgm'),('efeB','K:K16301','rgm'),('irtA','K:K14698','rgm'),('bfr','K:K03594','rgm'),('bfrB','N:bfrB','rgm'),
('znuA_FSM','K:K09815','fsm'),('zur','K:K09823','rgm'),('mntH','K:K03322','rgm'),('sitA','N:sitA','rgm'),('mbtB','K:K04788','rgm'),
('thiC','K:K03147','rgm'),('thiD','K:K00941','rgm'),('thiE','K:K00788','rgm'),('thiG','K:K03149','rgm'),('thiO','K:K03153','rgm'),('thiM_FSM','N:thiM','fsm'),('tenA','K:K20896','rgm'),
('panD','K:K01579','rgm'),('panE','K:K00077','rgm'),('panB','K:K00606','rgm'),('egtA','N:egtA','rgm'),('egtB','N:egtB','rgm'),('egtC','N:egtC','rgm'),('egtD','N:egtD','rgm'),('egtE','N:egtE','rgm'),
]
out=[]
for q,key,src in Q:
    order=pref if src=='rgm' else fsmpref+[g for g in members(s1,'FSM')]
    if src=='rgm': order=order+members(s1,'RGM_other')+members(s1,'RGM_sister')
    got=None
    for g in order:
        v=d['prot'].get((g,key))
        if v: got=(g,v[0][0]); break
    out.append((q,key,got[0] if got else 'NA',got[1] if got else 'NA'))
    print(q,key,got)
pd.DataFrame(out,columns=['query','key','genome','protein']).to_csv(F+'tblastn_queries.tsv',sep='\t',index=False)
