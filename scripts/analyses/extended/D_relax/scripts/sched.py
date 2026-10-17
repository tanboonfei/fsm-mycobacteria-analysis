# Minimal thread-budgeted job runner. jobs.tsv lines: threads<TAB>name<TAB>shell command (appendable while running).
# budget.txt holds the current thread budget. status in sched_status.tsv
import subprocess, time, os, sys
D=sys.argv[1]; jf=D+'/jobs.tsv'; bf=D+'/budget.txt'; sf=D+'/sched_status.tsv'
running={}; done=set()
def log(*a):
    with open(sf,'a') as o: o.write('\t'.join(map(str,a))+'\t'+time.strftime('%Y-%m-%d %H:%M:%S')+'\n')
while True:
    jobs=[l.rstrip('\n').split('\t',2) for l in open(jf) if l.strip() and not l.startswith('#')]
    try: budget=int(open(bf).read().strip())
    except Exception: budget=1
    for n,(p,t) in list(running.items()):
        r=p.poll()
        if r is not None: log('END',n,r); done.add(n); del running[n]
    used=sum(t for p,t in running.values())
    for th,n,cmd in jobs:
        th=int(th)
        if n in done or n in running: continue
        if used+th<=budget:
            p=subprocess.Popen(['bash','-c',cmd]); running[n]=(p,th); used+=th; log('START',n,th)
    if os.path.exists(D+'/STOP') and not running: break
    time.sleep(20)
