#!/usr/bin/env python3
"""Validate follow-up jobs and optionally submit each once."""
import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--campaign", choices=["all", "long_N4000", "N8000"], default="all")
    args = parser.parse_args()
    rows = json.loads((ROOT / ("smoke_jobs.json" if args.smoke else "jobs.json")).read_text())
    ledger = ROOT / ("smoke_submissions.jsonl" if args.smoke else "submissions.jsonl")
    history = [json.loads(line) for line in ledger.read_text().splitlines()] if ledger.exists() else []
    last = {record["directory"]: record for record in history}
    states = {}
    if args.retry_failed and last:
        ids = ",".join(record["job_id"] for record in last.values())
        output = subprocess.check_output(
            ["sacct", "-j", ids, "-n", "-P", "-X", "-o", "JobIDRaw,State"], text=True
        )
        states = dict(line.split("|")[:2] for line in output.splitlines() if line.strip())
    active = set(subprocess.check_output(
        ["squeue", "-h", "-u", os.environ["USER"], "-o", "%j"], text=True
    ).splitlines())
    pending = []
    for row in rows:
        if args.campaign != "all" and row["campaign"] != args.campaign:
            continue
        folder = Path(row["directory"])
        manifest = json.loads((folder / "manifest.json").read_text())
        assert manifest["campaign"] == row["campaign"]
        assert manifest["N"] == row["N"]
        for key in ["input_hashes", "data_hashes"]:
            for name, expected in manifest[key].items():
                assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == expected, folder / name
        if row["name"] in active:
            continue
        if str(folder) in last:
            if not args.retry_failed:
                continue
            state = (states.get(last[str(folder)]["job_id"], "").split() or [""])[0]
            # Retry only jobs that ended unsuccessfully; skip completed or unknown ones.
            if state not in {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL"}:
                continue
        subprocess.run(["bash", "-n", str(folder / "job.sbatch")], check=True)
        result = subprocess.run(
            ["sbatch", "--test-only", str(folder / "job.sbatch")],
            capture_output=True, text=True, check=False
        )
        assert result.returncode == 0, (folder, result.stderr)
        pending.append(row)
    print(f"Validated allocations: {len(pending)}", flush=True)
    if args.submit:
        with ledger.open("a") as output:
            for row in pending:
                job_id = subprocess.check_output(
                    ["sbatch", "--parsable", str(Path(row["directory"]) / "job.sbatch")],
                    text=True
                ).strip().split(";")[0]
                assert job_id.isdigit(), job_id
                output.write(json.dumps({**row, "job_id": job_id}) + "\n")
                output.flush()
                os.fsync(output.fileno())
                print("Submitted", job_id, row["name"], flush=True)


if __name__ == "__main__":
    main()
