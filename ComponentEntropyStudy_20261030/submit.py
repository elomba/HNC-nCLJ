#!/usr/bin/env python3
"""Validate each allocation then submit once, recording all actual job IDs."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--smoke',action='store_true');p.add_argument('--submit',action='store_true');p.add_argument('--retry-failed',action='store_true');a=p.parse_args()
rows=json.loads((ROOT/('smoke_jobs.json' if a.smoke else 'jobs.json')).read_text())
ledger=ROOT/('smoke_submissions.jsonl' if a.smoke else 'submissions.jsonl')
history=[json.loads(line) for line in ledger.read_text().splitlines()] if ledger.exists() else []
last={r['directory']:r for r in history}
known=set(last)
states={}
if a.retry_failed and last:
    ids=','.join(r['job_id'] for r in last.values())
    accounting=subprocess.check_output(['sacct','-j',ids,'-n','-P','-X','-o','JobIDRaw,State'],text=True)
    states=dict(line.split('|')[:2] for line in accounting.splitlines() if line.strip())
active=set(subprocess.check_output(['squeue','-h','-u',os.environ['USER'],'-o','%j'],text=True).splitlines())
pending=[]
for row in rows:
    d=Path(row['directory']);m=json.loads((d/'manifest.json').read_text())
    for group in ['input_hashes','data_hashes']:
        for path,digest in m[group].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
    if row['name'] in active:continue
    if str(d) in known:
        if not a.retry_failed:continue
        state=states.get(last[str(d)]['job_id'],'')
        assert state.split()[0] in ['FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL'],(d,state)
    subprocess.run(['bash','-n',str(d/'job.sbatch')],check=True)
    res=subprocess.run(['sbatch','--test-only',str(d/'job.sbatch')],capture_output=True,text=True)
    assert res.returncode==0,(d,res.stderr)
    pending.append(row)
print('Validated allocations:',len(pending),flush=True)
if a.submit:
    with ledger.open('a') as f:
        for row in pending:
            job_id=subprocess.check_output(['sbatch','--parsable',str(Path(row['directory'])/'job.sbatch')],text=True).strip().split(';')[0]
            assert job_id.isdigit(),job_id
            record=dict(row,job_id=job_id)
            f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
            print('Submitted',job_id,row['name'],flush=True)
