#!/bin/bash
#SBATCH --share
#SBATCH --partition=short
#SBATCH --job-name=goal_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#$ -N goal_pipeline


help_and_exit() {
  local retval=${1:-1}
  cat <<EOF
${0##*/} [-o|--output-dir] [-p|--profile] [-n|--no-resume] input-dir
EOF
  exit "$retval"
}

NO_RESUME="FALSE"
while (( $# )); do
  case $1 in
    -o|--output-dir)  OUTPUT_DIR=$2; shift ;;
    -p|--profile)     PROFILE=$2; shift ;;
    -n|--no-resume)     NO_RESUME="TRUE"; shift ;;
    -*)        printf 'Unknown option: %q\n\n' "$1"
               help_and_exit 1 ;;
    *)         args+=( "$1" ) ;;
  esac
  shift
done
set -- "${args[@]}"

INPUT_DIR=$1
INPUT_NAME="$(basename $INPUT_DIR)"
if [ ! -z ${OUTPUT_DIR+x} ];
then
  OUTPUT_DIR_OPT="--output $OUTPUT_DIR"
fi

if [ ! -z ${PROFILE+x} ];
then
  if [ "$PROFILE" == "uab_cheaha" ] ||  [ "$PROFILE" == "uab_cheaha_local" ]
  then
    module load Singularity
    module load Nextflow
    mkdir -p /scratch/$USER/GOAL/logs
    WORK_DIR_OPT="-w /scratch/$USER/GOAL/work_dir"
  fi
  if [ "$PROFILE" == "uab_hydrogen" ]
  then
    module load Nextflow
    WORK_DIR_OPT="-w /scratch/goal_pipeline/work_dir"
  fi
  if [ "$PROFILE" == "uab_local" ]
  then
    WORK_DIR_OPT="-w $PWD/work_dir"
  fi
  PROFILE_OPT="-profile $PROFILE"
else
  # Assume Cheaha
  module load Singularity
  module load Nextflow
  mkdir -p /scratch/$USER/GOAL/logs
  WORK_DIR_OPT="-w /scratch/$USER/GOAL/work_dir"
  PROFILE_OPT="-profile uab_cheaha"
fi

if [ "$NO_RESUME" == "FALSE" ]
then
  RESUME_OPT="-resume"
else
  RESUME_OPT=""
fi

nextflow run ./goalConsensus.nf \
  $WORK_DIR_OPT \
  --input $INPUT_DIR $OUTPUT_DIR_OPT \
  $PROFILE_OPT $RESUME_OPT \
  -with-report logs/${INPUT_NAME}/report.html \
  -with-timeline logs/${INPUT_NAME}/timeline.html \
  -with-trace logs/${INPUT_NAME}/trace.txt
