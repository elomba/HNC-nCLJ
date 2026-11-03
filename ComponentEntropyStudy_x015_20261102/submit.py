#!/usr/bin/env python3
"""Validate allocations locally; submit with --submit (add --smoke for smoke jobs); ledger-protected."""
import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def validate(smoke=False):
    rows = json.loads((ROOT / ('smoke_jobs.json' if smoke else 'jobs.json')).read_text())
    assert len(rows) == (8 if smoke else 32)
    for row in rows:
        folder = Path(row['directory'])
        manifest = json.loads((folder / 'manifest.json').read_text())
        assert len(manifest['densities']) == (2 if smoke else 19)
        assert manifest['seed'] in (41057, 52691, 67433, 83819)
        for group in ('input_hashes', 'data_hashes'):
            for path, expected in manifest[group].items():
                assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, path
        subprocess.run(['bash', '-n', str(folder / 'job.sbatch')], check=True)
        script = (folder / 'job.sbatch').read_text()
        assert 'dense-stage resume is unsupported' in script
        assert 'analyze_md.py' in script
    return rows


def slurm_test(rows):
    for row in rows:
        script = str(Path(row['directory']) / 'job.sbatch')
        result = subprocess.run(['sbatch', '--test-only', script], capture_output=True, text=True,
                                check=False)
        if result.returncode:
            raise RuntimeError(f'sbatch --test-only failed for {script}: {result.stderr.strip()}')
    print(f'Slurm test-only accepted: {len(rows)} allocations', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--submit', action='store_true')
    parser.add_argument('--slurm-test', action='store_true')
    args = parser.parse_args()
    rows = validate(args.smoke)
    print(f'Validated allocations: {len(rows)} (hashes, stage counts, shell syntax)', flush=True)
    if args.slurm_test or args.submit:
        slurm_test(rows)
    if args.submit:
        ledger = ROOT / ('smoke_submissions.jsonl' if args.smoke else 'submissions.jsonl')
        history = [json.loads(line) for line in ledger.read_text().splitlines()] if ledger.exists() else []
        submitted = {row['directory'] for row in history}
        with ledger.open('a') as output:
            for row in rows:
                if row['directory'] in submitted:
                    continue
                job_id = subprocess.check_output(
                    ['sbatch', '--parsable', str(Path(row['directory']) / 'job.sbatch')],
                    text=True).strip().split(';')[0]
                assert job_id.isdigit(), job_id
                output.write(json.dumps({**row, 'job_id': job_id}) + '\n')
                output.flush()
                os.fsync(output.fileno())
                print('Submitted', job_id, row['name'], flush=True)


if __name__ == '__main__':
    main()
