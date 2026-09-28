#!/usr/bin/env python3
"""
Generates one SLURM sbatch script per sample (fastp -> align -> filter ->
freebayes) and submits each with `sbatch`. Edit CONFIG, then:
    python submit_jobs.py
"""

import glob
import os
import subprocess

# ============================== CONFIG ==============================
RAWDATA_ROOT = "/hpcdata/Mimir/adrian/research/095/data/X204SC24122660-Z01-F001/01.RawData"
REFERENCE_FASTA = "/hpcdata/Mimir/adrian/research/095/ref/human_nac.nascent.fa"
JAVA_BIN = "/users/home/adrian/software/java/jdk-21.0.8/bin/java"
PICARD_JAR = "/users/home/adrian/software/picard/picard.jar"
FREEBAYES_BIN = "/users/home/adrian/software/freebayes/freebayes-1.3.10-linux-amd64-static"
CLEANED_DIR = "/hpcdata/Mimir/adrian/research/095/cleaned"
ALIGNED_DIR = "/hpcdata/Mimir/adrian/research/095/aligned"
CALLED_DIR = "/hpcdata/Mimir/adrian/research/095/called"
JOBSCRIPT_DIR = "/hpcdata/Mimir/adrian/research/095/jobscripts"

THREADS = 16
NTASKS_PER_NODE = 8

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

SAMPLE_ID="{sample_id}"
SAMPLE_DIR="{sample_dir}"
REFERENCE_FASTA="{reference_fasta}"
JAVA_BIN="{java_bin}"
PICARD_JAR="{picard_jar}"
FREEBAYES_BIN="{freebayes_bin}"
THREADS={threads}

mkdir -p "{cleaned_dir}" "{aligned_dir}" "{called_dir}"

# ---------- fastp: concatenate lanes, then trim ----------
CAT_R1="{cleaned_dir}/${{SAMPLE_ID}}_R1.cat.fastq.gz"
CAT_R2="{cleaned_dir}/${{SAMPLE_ID}}_R2.cat.fastq.gz"
TRIM_R1="{cleaned_dir}/${{SAMPLE_ID}}_R1.trim.fastq.gz"
TRIM_R2="{cleaned_dir}/${{SAMPLE_ID}}_R2.trim.fastq.gz"
FASTP_JSON="{cleaned_dir}/${{SAMPLE_ID}}_fastp.json"
FASTP_HTML="{cleaned_dir}/${{SAMPLE_ID}}_fastp.html"

mapfile -t R1_LANES < <(find "$SAMPLE_DIR" -name "*_1.fq.gz" | sort)
mapfile -t R2_LANES < <(find "$SAMPLE_DIR" -name "*_2.fq.gz" | sort)
echo "[$SAMPLE_ID] concatenating ${{#R1_LANES[@]}} lane(s)"
cat "${{R1_LANES[@]}}" > "$CAT_R1"
cat "${{R2_LANES[@]}}" > "$CAT_R2"

echo "[$SAMPLE_ID] fastp"
fastp -i "$CAT_R1" -I "$CAT_R2" -o "$TRIM_R1" -O "$TRIM_R2" \\
    --detect_adapter_for_pe \\
    --correction \\
    --cut_right --cut_mean_quality 25 \\
    --trim_poly_g \\
    --qualified_quality_phred 20 \\
    --length_required 36 \\
    --thread "$THREADS" \\
    -j "$FASTP_JSON" -h "$FASTP_HTML"

# ---------- align, calmd, dedup, filter ----------
SORTED_BAM="{aligned_dir}/${{SAMPLE_ID}}.sorted.bam"
CALMD_BAM="{aligned_dir}/${{SAMPLE_ID}}.calmd.bam"
DEDUP_BAM="{aligned_dir}/${{SAMPLE_ID}}.dedup.bam"
DUP_METRICS="{aligned_dir}/${{SAMPLE_ID}}.dup_metrics.txt"
FILTERED_BAM="{aligned_dir}/${{SAMPLE_ID}}.filtered.bam"
VCF_OUT="{called_dir}/${{SAMPLE_ID}}.vcf"

echo "[$SAMPLE_ID] aligning"
bwa mem -M -t "$THREADS" -R "@RG\\tID:${{SAMPLE_ID}}\\tSM:${{SAMPLE_ID}}\\tPL:ILLUMINA" \\
    "$REFERENCE_FASTA" "$TRIM_R1" "$TRIM_R2" \\
    | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
samtools index "$SORTED_BAM"

echo "[$SAMPLE_ID] calmd"
samtools calmd -b "$SORTED_BAM" "$REFERENCE_FASTA" > "$CALMD_BAM"
samtools index "$CALMD_BAM"

echo "[$SAMPLE_ID] MarkDuplicates"
"$JAVA_BIN" -Xmx24g -jar "$PICARD_JAR" MarkDuplicates \\
    I="$CALMD_BAM" O="$DEDUP_BAM" M="$DUP_METRICS"

echo "[$SAMPLE_ID] filtering"
samtools view -b -q 30 -f 2 -F 0xD00 "$DEDUP_BAM" > "$FILTERED_BAM"
samtools index "$FILTERED_BAM"

# ---------- freebayes, single-threaded ----------
echo "[$SAMPLE_ID] calling variants (freebayes, single-threaded)"
"$FREEBAYES_BIN" \\
    -f "$REFERENCE_FASTA" \\
    --pooled-continuous \\
    --min-alternate-count 1 \\
    --min-alternate-fraction 0 \\
    --no-population-priors \\
    --hwe-priors-off \\
    --allele-balance-priors-off \\
    -i -X -u \\
    --min-mapping-quality 30 \\
    --min-base-quality 20 \\
    --min-repeat-entropy 1 \\
    "$FILTERED_BAM" > "$VCF_OUT"

echo "[$SAMPLE_ID] done -> $VCF_OUT"
"""

os.makedirs(JOBSCRIPT_DIR, exist_ok=True)

sample_dirs = sorted(d for d in glob.glob(f"{RAWDATA_ROOT}/*") if os.path.isdir(d))
print(f"Found {len(sample_dirs)} sample directories")
if not sample_dirs:
    raise SystemExit(f"ERROR: no sample directories found under RAWDATA_ROOT={RAWDATA_ROOT}")

for sample_dir in sample_dirs:
    sample_id = os.path.basename(sample_dir)

    script_content = JOB_TEMPLATE.format(
        sample_id=sample_id,
        sample_dir=sample_dir,
        reference_fasta=REFERENCE_FASTA,
        java_bin=JAVA_BIN,
        picard_jar=PICARD_JAR,
        freebayes_bin=FREEBAYES_BIN,
        cleaned_dir=CLEANED_DIR,
        aligned_dir=ALIGNED_DIR,
        called_dir=CALLED_DIR,
        jobscript_dir=JOBSCRIPT_DIR,
        threads=THREADS,
        ntasks_per_node=NTASKS_PER_NODE,
    )

    script_path = os.path.join(JOBSCRIPT_DIR, f"{sample_id}.sh")
    with open(script_path, "w") as fh:
        fh.write(script_content)
    print(f"Wrote {script_path}")

    if SUBMIT_JOBS:
        result = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
        print("  submitted:", result.stdout.strip() or result.stderr.strip())

if SUBMIT_JOBS:
    print(f"\nSubmitted {len(sample_dirs)} jobs. Check with: squeue -u $USER")
else:
    print(f"\nWrote {len(sample_dirs)} job scripts to {JOBSCRIPT_DIR} -- nothing submitted yet.")
    print("Inspect them (e.g. `cat {}/PN1_Rep1.sh`), then set SUBMIT_JOBS = True and rerun to submit all.".format(JOBSCRIPT_DIR))