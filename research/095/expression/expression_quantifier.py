#!/usr/bin/env python3
"""
Generate one SLURM sbatch script per sample and submit them for kallisto quant.
No fastp step here -- reuses already-trimmed FASTQs from the earlier fastp run.

Follows the same write-then-submit pattern as the project's other job generator:
set SUBMIT_JOBS = False first, inspect the generated .sbatch files, then set
SUBMIT_JOBS = True and rerun to actually submit.
"""

import glob
import os
import sys
import subprocess

# ==================== EDIT THIS SECTION ====================

TRIMMED_FASTQ_DIR = "/hpcdata/Mimir/adrian/research/095/cleaned"
R1_PATTERN = "*_R1.trim.fastq.gz"
R1_SUFFIX = "_R1"     # substring that separates sample_id from the R1/R2 marker
R2_SUFFIX = "_R2"

KALLISTO_INDEX = "/hpcdata/Mimir/adrian/research/095/ref/index.idx"
OUT_DIR = "/hpcdata/Mimir/adrian/research/095/quant"
JOBSCRIPT_DIR = "/hpcdata/Mimir/adrian/research/095/jobscripts_quant"

EXPECTED_SAMPLE_COUNT = 27   # sanity check only; script still runs if this doesn't match, just warns

THREADS = 16          # threads handed to kallisto's -t flag
NTASKS_PER_NODE = 8   # matches the existing project convention for this cluster

# Any module-load or environment setup your cluster needs before "kallisto" is on PATH.
# Leave as empty string if kallisto is already available in your login/job environment.
EXTRA_SETUP_CMDS = ""   # e.g. "module load kallisto/0.50.1"

KALLISTO_BIN = "kallisto"

SUBMIT_JOBS = True   # first pass: False (just writes the .sh files for you to inspect).
                       # Once you've checked them, set this to True and rerun to submit.
# ======================================================================

JOB_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={sample_id}
#SBATCH --partition=mimir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --hint=nomultithread
#SBATCH --output={jobscript_dir}/{sample_id}.log
#SBATCH --error={jobscript_dir}/{sample_id}.err

set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

{extra_setup}

SAMPLE_ID="{sample_id}"
KALLISTO_INDEX="{index}"
THREADS={threads}

mkdir -p "{out_dir}/{sample_id}"

{kallisto_bin} quant \\
    -i "$KALLISTO_INDEX" \\
    -o "{out_dir}/{sample_id}" \\
    -t $THREADS \\
    "{r1}" "{r2}"
"""


def print_config():
    print("=" * 70)
    print("CONFIG")
    print("=" * 70)
    print(f"TRIMMED_FASTQ_DIR   = {TRIMMED_FASTQ_DIR}")
    print(f"R1_PATTERN          = {R1_PATTERN}")
    print(f"R1_SUFFIX / R2_SUFFIX = {R1_SUFFIX} / {R2_SUFFIX}")
    print(f"KALLISTO_INDEX      = {KALLISTO_INDEX}")
    print(f"OUT_DIR             = {OUT_DIR}")
    print(f"JOBSCRIPT_DIR       = {JOBSCRIPT_DIR}")
    print(f"THREADS             = {THREADS}")
    print(f"NTASKS_PER_NODE     = {NTASKS_PER_NODE}")
    print(f"EXTRA_SETUP_CMDS    = {EXTRA_SETUP_CMDS!r}")
    print(f"SUBMIT_JOBS         = {SUBMIT_JOBS}")
    print("=" * 70)


def find_samples():
    pattern = os.path.join(TRIMMED_FASTQ_DIR, R1_PATTERN)
    print(f"\nGlobbing: {pattern}")
    r1_files = sorted(glob.glob(pattern))
    if not r1_files:
        sys.exit(f"No R1 files matched {pattern} -- check TRIMMED_FASTQ_DIR / R1_PATTERN.")
    print(f"Matched {len(r1_files)} R1 files.\n")

    samples = []
    for r1 in r1_files:
        r2 = r1.replace(R1_SUFFIX, R2_SUFFIX)
        sample_id = os.path.basename(r1).split(R1_SUFFIX)[0]
        if not os.path.exists(r2):
            print(f"  [{sample_id}] R1: {r1}")
            print(f"  [{sample_id}] R2: {r2}  <-- NOT FOUND, skipping")
            print()
            continue
        print(f"  [{sample_id}] R1: {r1}")
        print(f"  [{sample_id}] R2: {r2}")
        print()
        samples.append((sample_id, r1, r2))
    return samples


def write_job_scripts(samples):
    os.makedirs(JOBSCRIPT_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    print("\n" + "=" * 70)
    print("GENERATING JOB SCRIPTS")
    print("=" * 70)

    script_paths = []
    for sample_id, r1, r2 in samples:
        content = JOB_TEMPLATE.format(
            sample_id=sample_id,
            jobscript_dir=JOBSCRIPT_DIR,
            threads=THREADS,
            ntasks_per_node=NTASKS_PER_NODE,
            extra_setup=EXTRA_SETUP_CMDS,
            out_dir=OUT_DIR,
            kallisto_bin=KALLISTO_BIN,
            index=KALLISTO_INDEX,
            r1=r1,
            r2=r2,
        )
        script_path = os.path.join(JOBSCRIPT_DIR, f"{sample_id}.sh")

        kallisto_cmd = (
            f"{KALLISTO_BIN} quant -i {KALLISTO_INDEX} "
            f"-o {OUT_DIR}/{sample_id} -t {THREADS} {r1} {r2}"
        )
        print(f"\n[{sample_id}] writing {script_path}")
        print(f"[{sample_id}] kallisto command: {kallisto_cmd}")
        print(f"[{sample_id}] full job script:")
        print("-" * 70)
        print(content)
        print("-" * 70)

        with open(script_path, "w") as fh:
            fh.write(content)
        os.chmod(script_path, 0o755)
        script_paths.append((sample_id, script_path))
    return script_paths


def submit_jobs(script_paths):
    print("\n" + "=" * 70)
    print("SUBMITTING" if SUBMIT_JOBS else "DRY RUN (SUBMIT_JOBS = False, nothing will be submitted)")
    print("=" * 70)
    for sample_id, script_path in script_paths:
        if not SUBMIT_JOBS:
            print(f"[dry-run] would run: sbatch {script_path}")
            continue
        cmd = ["sbatch", script_path]
        print(f"[{sample_id}] running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[{sample_id}] FAILED to submit: {result.stderr.strip()}", file=sys.stderr)
        else:
            print(f"[{sample_id}] {result.stdout.strip()}")


def main():
    print_config()
    samples = find_samples()
    print(f"Found {len(samples)} samples with matched R1/R2 pairs.")
    if len(samples) != EXPECTED_SAMPLE_COUNT:
        print(f"WARNING: expected {EXPECTED_SAMPLE_COUNT} samples, found {len(samples)}. "
              f"Check TRIMMED_FASTQ_DIR / R1_PATTERN before trusting this batch.", file=sys.stderr)

    script_paths = write_job_scripts(samples)
    print(f"\nWrote {len(script_paths)} job scripts to {JOBSCRIPT_DIR}")

    submit_jobs(script_paths)


if __name__ == "__main__":
    main()