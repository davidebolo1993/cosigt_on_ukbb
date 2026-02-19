# COSIGT on UK Biobank Research Analysis Platform

This guide describes how to run the main steps of the [cosigt](https://github.com/davidebolo1993/cosigt) genotyping pipeline on the [UK Biobank Research Analysis Platform (RAP)](https://ukbiobank.dnanexus.com/) using the DNAnexus `app-swiss-army-knife` and Docker containers.

## Overview

The pipeline runs in three sequential stages:

| Step | Script | Runs | Instance |
|------|--------|------|----------|
| 1. Preprocess reference | `submit-rap-reference.sh` | Once | `mem3_ssd2_v2_x8` |
| 2. Preprocess samples | `submit-rap-sample.sh` | One job per sample (parallel) | `mem1_ssd1_v2_x8` |
| 3. Process regions | `submit-rap-region.sh` | One job per region (parallel) | `mem2_ssd1_v2_x16` |

Large input files (reference FASTA, CRAMs/BAMs, graph `.og` files) are accessed directly via the `/mnt/project` FUSE mount

## Prerequisites

- DNAnexus CLI installed and authenticated:

```bash
pip install dxpy
dx login
dx select <project>
dx cd <User>
```

- The following files uploaded (`dx upload`) in your **local working directory** before submitting:

```text
config.yaml
aln_map.tsv
graph_map.tsv
ref.fa #(or ref.fa.gz)
ref.fa.fai # (or ref.fa.gz.fai and ref.fa.gz.gzi)
chrA_START1_END1.og
chrA_START2_END2.og
chrN_STARTM_ENDM.og
SAMPLE1.bam #(or SAMPLE1.cram)
SAMPLE1.bam.bai #(or SAMPLE1.cram.crai)
SAMPLE2.bam #(or SAMPLE2.cram)
SAMPLE2.bam.bai #(or SAMPLE2.cram.crai)
SAMPLEN.bam #(or SAMPLEN.cram)
SAMPLEN.bam.bai #(or SAMPLEN.cram.crai)
```

## Configuration

### config.yaml

```yaml
input:         /mnt/project/<User>/cosigt_results           # where to READ pre-existing results from (unmapped reads, meryl database)
output:        ./cosigt_results                              # where to WRITE new outputs (local to the job, then uploaded)
reference:     /mnt/project/<User>/ref.fa
alignment_map: /mnt/project/<User>/aln_map.tsv
graph_map:     /mnt/project/<User>/graph_map.tsv

samples:
- SAMPLE1
- SAMPLE2
- SAMPLEN

regions:
- chrA_START1_END1
- chrA_START2_END2
- chrN_STARTM_ENDM
```

- `input` points to the **read-only** `/mnt/project` location where results from Steps 1 and 2 are generated. This is necessary for Step 3, since this reads the meryl DB and unmapped FASTAs without being able to write to `/mnt/project`.
- `output` is a **local path** relative to the job working directory. DNAnexus uploads everything written here after the job completes.

### aln_map.tsv

Tab-separated, alignment path first and sample name second (or viceversa)

```tsv
/mnt/project/<User>/SAMPLE1.bam	SAMPLE1
/mnt/project/<User>/SAMPLE2.bam	SAMPLE2
/mnt/project/<User>/SAMPLEN.bam	SAMPLEN
```

BAM/CRAM index files (`.bai`, `.csi`, `.crai`) must exist alongside the alignment files in the project.

### graph_map.tsv

Tab-separated, graph path first, region name second (or viceversa)

```tsv
/mnt/project/<User>/chrA_START1_END1.og chrA_START1_END1
/mnt/project/<User>/chrA_START2_END2.og chrA_START2_END2
/mnt/project/<User>/chrN_STARTM_ENDM.og chrN_STARTM_ENDM
```

## Step 0 — Pre-build all the docker images

```bash
docker pull davidebolo1993/cosigt-preprocess-reference:latest
docker save davidebolo1993/cosigt-preprocess-reference:latest | gzip > cosigt-preprocess-reference.tar.gz

docker pull davidebolo1993/cosigt-preprocess-sample:latest
docker save davidebolo1993/cosigt-preprocess-sample:latest | gzip > cosigt-preprocess-sample.tar.gz

docker pull davidebolo1993/cosigt-process-region:latest
docker save davidebolo1993/cosigt-process-region:latest | gzip > cosigt-process-region.tar.gz
```

These files are passed via `-iin` and loaded inside the job with `docker load` before `docker run`.

## Step 1 — Preprocess reference (once per reference genome)

```bash
dx run app-swiss-army-knife \
  -iin="config.yaml" \
  -iin="cosigt-preprocess-reference.tar.gz" \
  --ignore-reuse \
  -icmd="
  docker load < cosigt-preprocess-reference.tar.gz && \
  docker run --rm \
  -v ./:/work \
  -v /mnt/project:/mnt/project \
  davidebolo1993/cosigt-preprocess-reference:latest \
  config.yaml \
  8 \
  "     
  --destination="./"
  --name="preprocess-reference" --instance-type="mem3_ssd2_v2_x8"
```

**Output:** a Meryl k-mer database written to `./cosigt_results/meryl/reference/`, which DNAnexus uploads to your project destination.

## Step 2 — Preprocess samples (one job per sample)

The script reads the `samples:` block from `config.yaml` and submits **one independent job per sample**.

```bash
JOB_IDS=()
SAMPLES=$(awk '/^samples:/{flag=1; next} /^[a-zA-Z]/{flag=0} flag && /^- /{gsub(/^- /, ""); print}' config.yaml)

for SAMPLE in $SAMPLES; do
    echo "Submitting: $SAMPLE"
    JOB_ID=$(dx run app-swiss-army-knife \
    -iin="config.yaml" \
    -iin="aln_map.tsv" \
    -iin="cosigt-preprocess-sample.tar.gz" \
    --ignore-reuse \
    -icmd="
    docker load < cosigt-preprocess-sample.tar.gz && \
    docker run --rm \
    -v ./:/work \
    -v /mnt/project:/mnt/project \
    davidebolo1993/cosigt-preprocess-sample:latest \
    ${SAMPLE} \
    config.yaml \
    8" \
    --destination="./" \
    --instance-type="mem1_ssd1_v2_x8" \
    --name="preprocess-${SAMPLE}" \
    --brief \
    --yes)
    JOB_IDS+=("$JOB_ID")
done
echo "Job IDs:"
printf '%s\n' "${JOB_IDS[@]}"
```

**Output:** `./cosigt_test/samtools/fasta/${SAMPLE}/unmapped.fasta.gz` per sample, uploaded to your project.

Job IDs are printed at the end so one can monitor progress:

```bash
dx watch <job-id>
```

## Step 3 — Process regions (one job per region)

**Wait for Steps 1 and 2 to complete before running this step.**

Once all preprocessing outputs are uploaded one can run the final step:

```bash
JOB_IDS=()
REGIONS=$(awk '/^regions:/{flag=1; next} /^[a-zA-Z]/{flag=0} flag && /^- /{gsub(/^- /, ""); print}' config.yaml)

for REGION in $REGIONS; do
    echo "Submitting: $REGION"
    JOB_ID=$(dx run app-swiss-army-knife \
    -iin="config.yaml" \
    -iin="aln_map.tsv" \
    -iin="graph_map.tsv" \
    -iin="cosigt-process-region.tar.gz" \
    --ignore-reuse \
    -icmd="
    docker load < cosigt-process-region.tar.gz && \
    docker run --rm \
    -v ./:/work \
    -v /mnt/project:/mnt/project \
    davidebolo1993/cosigt-process-region:latest \
    ${REGION} \
    config.yaml \
    12" \
    --destination="./" \
    --instance-type="mem2_ssd1_v2_x16" \
    --name="process-${REGION}" \
    --brief \
    --yes)
    JOB_IDS+=("$JOB_ID")
done
echo "Job IDs:"
printf '%s\n' "${JOB_IDS[@]}"
```

**Output:** `./cosigt_test/cosigt/${SAMPLE}/${CHROM}/${REGION}/${REGION}.cosigt_genotype.tsv` per sample/region, uploaded to your project.

## Monitoring jobs

```bash
# List all running jobs
dx find jobs --state running

# Stream log for a specific job
dx watch <job-id>

# Wait until a job finishes
dx wait <job-id>
```

## Output structure

After all three steps complete, your DNAnexus project destination will contain:

```text
cosigt_test/
├── meryl/
│   └── reference/                  # Meryl k-mer DB (Step 1)
├── samtools/fasta/
│   ├── HG00096/
│   │   └── unmapped.fasta.gz       # Unmapped reads (Step 2)
│   └── HG00171/
│       └── unmapped.fasta.gz
└── cosigt/
    ├── HG00096/
    │   └── chr1/
    │       └── chr1_103304997_103901127/
    │           └── chr1_103304997_103901127.cosigt_genotype.tsv
    └── HG00171/
        └── chr1/
            └── chr1_103304997_103901127/
                └── chr1_103304997_103901127.cosigt_genotype.tsv
```
