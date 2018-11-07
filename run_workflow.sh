#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_workflow.sh --samples FILE --genome-dir DIR --ref-flat FILE \
  --jum-dir DIR --design FILE --output-dir DIR [options]

Required sample-sheet columns: sample, condition, read1, read2 (tab-separated).
The sheet must contain exactly two conditions.

Options:
  --threads N          STAR and JUM threads (default: 8)
  --read-length N      Read length supplied to JUM (default: 100)
  --cutoff VALUE       JUM p-value cutoff (default: 0.05)
  --read-files-command CMD  STAR decompression command, for example zcat
  --execute            Run commands; otherwise print a dry-run plan
  -h, --help           Show this help
EOF
}

samples_file=""
genome_dir=""
ref_flat=""
jum_dir=""
design_file=""
output_dir=""
threads=8
read_length=100
cutoff=0.05
read_files_command=""
execute=false

while (($#)); do
  case "$1" in
    --samples) samples_file=${2:?missing value}; shift 2 ;;
    --genome-dir) genome_dir=${2:?missing value}; shift 2 ;;
    --ref-flat) ref_flat=${2:?missing value}; shift 2 ;;
    --jum-dir) jum_dir=${2:?missing value}; shift 2 ;;
    --design) design_file=${2:?missing value}; shift 2 ;;
    --output-dir) output_dir=${2:?missing value}; shift 2 ;;
    --threads) threads=${2:?missing value}; shift 2 ;;
    --read-length) read_length=${2:?missing value}; shift 2 ;;
    --cutoff) cutoff=${2:?missing value}; shift 2 ;;
    --read-files-command) read_files_command=${2:?missing value}; shift 2 ;;
    --execute) execute=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'error: unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

for value_name in samples_file genome_dir ref_flat jum_dir design_file output_dir; do
  if [[ -z ${!value_name} ]]; then
    printf 'error: --%s is required\n' "${value_name//_/-}" >&2
    exit 2
  fi
done
[[ $threads =~ ^[1-9][0-9]*$ ]] || { echo 'error: --threads must be positive' >&2; exit 2; }
[[ $read_length =~ ^[1-9][0-9]*$ ]] || { echo 'error: --read-length must be positive' >&2; exit 2; }
[[ -f $samples_file ]] || { echo "error: sample sheet not found: $samples_file" >&2; exit 2; }
[[ -d $genome_dir ]] || { echo "error: genome directory not found: $genome_dir" >&2; exit 2; }
[[ -f $ref_flat ]] || { echo "error: refFlat file not found: $ref_flat" >&2; exit 2; }
[[ -d $jum_dir ]] || { echo "error: JUM directory not found: $jum_dir" >&2; exit 2; }
[[ -f $design_file ]] || { echo "error: design file not found: $design_file" >&2; exit 2; }

for required in JUM_A.sh JUM_B.sh JUM_C.sh R_script_JUM.R; do
  [[ -f $jum_dir/$required ]] || { echo "error: missing JUM file: $jum_dir/$required" >&2; exit 2; }
done

samples=()
conditions=()
read1_files=()
read2_files=()
declare -A seen_samples=()
declare -A condition_counts=()

IFS= read -r header < "$samples_file"
[[ $header == $'sample\tcondition\tread1\tread2' ]] || {
  echo 'error: sample sheet header must be: sample<TAB>condition<TAB>read1<TAB>read2' >&2
  exit 2
}

while IFS=$'\t' read -r sample condition read1 read2 extra; do
  [[ -n $sample || -n $condition || -n $read1 || -n $read2 ]] || continue
  [[ -z ${extra:-} && -n $sample && -n $condition && -n $read1 && -n $read2 ]] || {
    echo 'error: each sample row must have exactly four non-empty fields' >&2
    exit 2
  }
  [[ $sample =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "error: unsafe sample name: $sample" >&2; exit 2; }
  [[ $condition =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "error: unsafe condition name: $condition" >&2; exit 2; }
  [[ -z ${seen_samples[$sample]+x} ]] || { echo "error: duplicate sample: $sample" >&2; exit 2; }
  [[ -f $read1 && -f $read2 ]] || { echo "error: reads not found for sample: $sample" >&2; exit 2; }
  seen_samples[$sample]=1
  condition_counts[$condition]=$(( ${condition_counts[$condition]:-0} + 1 ))
  samples+=("$sample")
  conditions+=("$condition")
  read1_files+=("$read1")
  read2_files+=("$read2")
done < <(tail -n +2 "$samples_file")

((${#samples[@]} > 1)) || { echo 'error: at least two samples are required' >&2; exit 2; }
((${#condition_counts[@]} == 2)) || { echo 'error: exactly two conditions are required' >&2; exit 2; }
mapfile -t condition_names < <(printf '%s\n' "${!condition_counts[@]}" | LC_ALL=C sort)
condition1=${condition_names[0]}
condition2=${condition_names[1]}

run_cmd() {
  printf '+ '
  printf '%q ' "$@"
  printf '\n'
  if $execute; then "$@"; fi
}

run_in_dir() {
  local directory=$1
  shift
  printf '+ (cd %q && ' "$directory"
  printf '%q ' "$@"
  printf ')\n'
  if $execute; then (cd "$directory" && "$@"); fi
}

join_by_comma() {
  local first=true value
  for value in "$@"; do
    if $first; then first=false; else printf ','; fi
    printf '%s_' "$value"
  done
}

read_command_args=()
if [[ -n $read_files_command ]]; then
  read_command_args=(--readFilesCommand "$read_files_command")
fi

pass1=$output_dir/pass1
pass2=$output_dir/pass2
analysis=$output_dir/analysis
run_cmd mkdir -p "$pass1" "$pass2" "$analysis"

junctions=()
condition1_samples=()
condition2_samples=()
for index in "${!samples[@]}"; do
  sample=${samples[$index]}
  condition=${conditions[$index]}
  prefix1=$pass1/${sample}_
  prefix2=$pass2/${sample}_
  run_cmd STAR --runThreadN "$threads" --genomeDir "$genome_dir" \
    --outFileNamePrefix "$prefix1" --readFilesIn "${read1_files[$index]}" "${read2_files[$index]}" \
    --outSJfilterReads Unique "${read_command_args[@]}"
  junctions+=("${prefix1}SJ.out.tab")
  if [[ $condition == "$condition1" ]]; then
    condition1_samples+=("$sample")
  else
    condition2_samples+=("$sample")
  fi
done

for index in "${!samples[@]}"; do
  sample=${samples[$index]}
  prefix2=$pass2/${sample}_
  run_cmd STAR --runThreadN "$threads" --genomeDir "$genome_dir" \
    --outFileNamePrefix "$prefix2" --readFilesIn "${read1_files[$index]}" "${read2_files[$index]}" \
    --outSJfilterReads Unique "${read_command_args[@]}" --outSAMstrandField intronMotif --outFilterMultimapNmax 1 \
    --sjdbFileChrStartEnd "${junctions[@]}" --outSAMtype BAM Unsorted
  run_cmd samtools sort -@ "$threads" -o "$analysis/${sample}_Aligned.out_sorted.bam" \
    "${prefix2}Aligned.out.bam"
  run_cmd samtools index "$analysis/${sample}_Aligned.out_sorted.bam"
done

condition1_arg=$(join_by_comma "${condition1_samples[@]}")
condition2_arg=$(join_by_comma "${condition2_samples[@]}")
count1=${condition_counts[$condition1]}
count2=${condition_counts[$condition2]}
threshold1=$((count1 > 1 ? count1 - 1 : 1))
threshold2=$((count2 > 1 ? count2 - 1 : 1))
total=$((count1 + count2))

run_in_dir "$analysis" bash "$jum_dir/JUM_A.sh" --Folder "$jum_dir" \
  --JuncThreshold 5 --Condition1_fileNum_threshold "$threshold1" \
  --Condition2_fileNum_threshold "$threshold2" --IRthreshold 5 \
  --Readlength "$read_length" --Thread "$threads" \
  --Condition1SampleName "$condition1_arg" --Condition2SampleName "$condition2_arg"

jum_diff=$analysis/JUM_diff
run_in_dir "$jum_diff" Rscript "$jum_dir/R_script_JUM.R" "$design_file"
run_in_dir "$jum_diff" bash "$jum_dir/JUM_B.sh" --Folder "$jum_dir" \
  --Test pvalue --Cutoff "$cutoff" --TotalFileNum "$total" \
  --Condition1_fileNum_threshold "$threshold1" --Condition2_fileNum_threshold "$threshold2" \
  --Condition1SampleName "$condition1_arg" --Condition2SampleName "$condition2_arg"

final_dir=$jum_diff/FINAL_JUM_OUTPUT_pvalue_${cutoff}
run_in_dir "$final_dir" bash "$jum_dir/JUM_C.sh" --Folder "$jum_dir" \
  --Test pvalue --Cutoff "$cutoff" --TotalCondition1FileNum "$count1" \
  --TotalCondition2FileNum "$count2" --REF "$ref_flat"

if ! $execute; then
  echo 'Dry run only. Re-run with --execute after reviewing the plan.'
fi
