"""Habitat summary of FSM-specific 16S records (SILVA 144 SSURef), excluding four FSM-MAG 16S that are not mycobacterial."""
import csv,re,collections
bad={'spire_mag_01108461_16S1','spire_mag_01753944_16S1','spire_mag_01960731_16S1','spire_mag_02102612_16S1'}
R=[r for r in csv.DictReader(open('silva16S_hits_annotated.tsv'),delimiter='\t') if r['query_is_FSM']=='True' and r['FSM_specific']=='True' and r['best_query'] not in bad]
def cat(r):
    t=' '.join(r[k] for k in ('isolation_source','host','definition','title')).lower()
    if re.search(r'amoeb|protist|protozo|ciliat|endosymb|symbiont|co-?cultur',t): return 'protist/symbiont context'
    if 'homo sapiens' in t or 'skin' in t: return 'human skin'
    if re.search(r'wastewater|sludge',t): return 'wastewater/activated sludge'
    if re.search(r'sediment',t): return 'sediment'
    if re.search(r'moss',t): return 'aquatic moss (lake)'
    if re.search(r'lake|reservoir|river|freshwater|epilimnion|water',t): return 'freshwater lake/reservoir/river water'
    if re.search(r'caspian',t): return 'brackish water column'
    return 'unspecified metagenome'
seen={}
for r in R: seen.setdefault(r['accession'],r)
c=collections.Counter(cat(r) for r in seen.values())
with open('fsm16S_specific_records.tsv','w') as o:
    o.write('accession\tbest_FSM_query\tpct_identity\tpanel_best_FSM\tpanel_best_nonFSM\tcategory\tisolation_source\thost\tcountry\torganism_label\ttitle\n')
    for a,r in sorted(seen.items(),key=lambda x:cat(x[1])):
        o.write('\t'.join([a,r['best_query'],r['best_pid'],r['panel_best_FSM'],r['panel_best_nonFSM'],cat(r),r['isolation_source'],r['host'],r['country'],r['organism'],r['title']])+'\n')
print(len(R),'hits;',len(seen),'records;',sum(float(r['best_pid'])>=99 for r in seen.values()),'at >=99%')
for k,v in c.most_common(): print(v,k)
