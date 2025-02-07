# Software for Clinical Health Omics Oncology Laboratories
School is a collection of genomics analysis workflows that are used for detecting single nucleotide variants (SNVs), insertions/deletions (indels), copy number variants (CNVs) and translocations from RNA and DNA sequencing. These workflows were first developed and validated in a CLIA laboratory at UTSW, and will continue to be developed and maintained by the Genomics Organization for Academic Laboratories (GOAL) Consortium.

## Prerequisites

These bioinformatics pipelines use **`nextflow`**, a framework for defining and executing a directed acyclic graph (DAG) of interdependent steps. Also required is either **`singularity`** or **`docker`**, executors for tools that are [containerized](https://www.docker.com/resources/what-container) for portability across computing environments.

[Follow these instructions](https://gist.github.com/ckandoth/982ce140b4dd9d6bf72a780c05a549a3) to install Nextflow and Singularity in a Linux environment. For better portability across computing environments (Linux, macOS, Windows), [follow these instructions](https://docs.docker.com/get-docker) to install Docker. Docker requires administrative rights, which you normally have on a personal laptop/workstation. But in shared computers like HPC clusters, there are [valid concerns](https://duo.com/decipher/docker-bug-allows-root-access-to-host-file-system) against installing Docker, and then Singularity makes more sense.

Some HPC clusters will have these tools pre-installed as [environment modules](https://modules.readthedocs.io/en/latest/). Use command `module avail` to see what's available, and `module load` to load them into your `$PATH`. But make sure you have Nextflow 20.07.1 or newer and Singularity 3.5.2 or newer.

## Quick Start

Clone a branch of this repo that you want to test and `cd` into it:
```bash
git clone -b UFHPL_branch_1 --single-branch https://github.com/goalconsortium/goal_school.git
cd goal_school
```

Download and unzip resource files needed by the pipeline (39GB download that unzips to 43GB):
```bash
curl -LO https://reference-files-bucket.s3.amazonaws.com/Reference_Files.zip
unzip Reference_Files.zip
```

In a folder named `fastq`, download small FASTQs created from DNA-seq of a tumor (or use your own FASTQs):
```bash
mkdir fastq
wget -P fastq https://github.com/mskcc/roslin-variant/raw/2.4.x/setup/examples/data/fastq/DU874145-T/DU874145-T_IGO_00000_TEST_L001_R{1,2}_001.fastq.gz
```

Prepare a design file for the DNA-seq variant calling pipeline:
```bash
echo -e "SampleID\tCaseID\tTumorID\tNormalID\tFqR1\tFqR2" > fastq/design.txt
echo -e "DU874145-T\tDU874145\tDU874145-T\t\tDU874145-T_IGO_00000_TEST_L001_R1_001.fastq.gz\tDU874145-T_IGO_00000_TEST_L001_R2_001.fastq.gz" >> fastq/design.txt
```

Run the `goalConsensus.nf` pipeline using the `standard` profile that uses Nextflow with singularity on the local machine:
```bash
nextflow run -work-dir .nextflow_work -profile standard goalConsensus.nf --input fastq --output analysis --repoDir ${PWD} --seqrunid H7YRLADXX --genome Reference_Files
```

To run this workflow on a different computing environment, lookup the institute-specific profiles in `nextflow.config` and/or create your own.

# Run Nextflow Workflows

This workflow can run either DNA or RNA sequencing. Please determine the desired configuration to achieve the proper analysis run.

## SlideRule: DNA Workflow

### DNA Design File

The design file must named design.txt and be in tab seperated format for the workflows. This workflow can be run with tumor-only or with tumor and normal pairs. If running tumor-only then do not include the NormalID collumn in the design file.

| SampleID | CaseID | TumorID | NormalID | FqR1 | FqR2 |
|---|---|---|---|---|---|
| Sample1 | Fam1 | Sample1 | Sample2 | Sample1.R1.fastq.gz | Sample1.R2.fastq.gz |
| Sample2 | Fam1 | Sample1 | Sample2 | Sample2.R1.fastq.gz | Sample2.R2.fastq.gz |
| Sample3 | Fam2 | Sample3 | Sample4 | Sample3.R1.fastq.gz | Sample3.R2.fastq.gz |
| Sample4 | Fam2 | Sample3 | Sample4 | Sample4.R1.fastq.gz | Sample4.R2.fastq.gz |

### DNA Parameters
* **--input**
  * directory containing the design file and fastq files
  * default is set to *'${basedir}/fastq'*
  * eg: **--input '/project/shared/bicf_workflow_ref/workflow_testdata/germline_variants/fastq'**
* **--output**
  * directory for the analysis output
  * default is set to *'${basedir}/analysis'*
  * eg: **--output '${basedir}/output'**
* **--pon**
  * pon file for mutect capture kit bed file
  * default is set to not included. If a pon file is included, then mutect will run with reference to the inputted file
  * eg: **--pon '/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/mutect2.pon.vcf.gz'** 
* **--capturedir**
  * directory containing capture bed and supporting files
  * eg: **--capturedir '/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme'**
* **--markdups**
  * select the method for marking duplicates from *'picard'*, *'samtools'*, *'fgbio_umi'*, *'picard_umi'*, or *'none'*
  * default is set to *'fgbio_umi'*
  * eg: **--markdups 'picard_umi'**
* **--snpeff_vers**
  * version of reference genome for snpeff tool
  * default is set to *'GRCh38.86'*
  * eg: **--snpeff_vers 'GRCh38.86'**

### DNA Run Workflow

On Cheaha
```
RUN=<Sequencer Run ID>
sbatch analyze.job /data/user/$USER/GOAL/PIPELINE_INPUT/${RUN} -o /data/user/$USER/GOAL/PIPELINE_OUTPUT/${RUN}
```

On Hydrogen
```
RUN=<Sequencer Run ID>
qsub -cwd -N goal_pipeline analyze.job /scratch/goal_pipeline/PIPELINE_INPUT/${RUN} -o /scratch/PIPELINE_OUTPUT/GOAL/${RUN} -p uab_hydrogen
```

On local machine
```
RUN=<Sequencer Run ID>
./analyze.job <Path to input>/${RUN} -o <Path to output>/${RUN} -p uab_local
```
