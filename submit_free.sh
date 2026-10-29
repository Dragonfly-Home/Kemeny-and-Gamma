#!/bin/bash
#SBATCH --job-name=kemeny_free
#SBATCH --account=free
#SBATCH --partition=freecpuQ
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=kemeny-free-%j.out
#SBATCH --error=kemeny-free-%j.err

set -eo pipefail
trap 'rc=$?; printf "ERROR: line=%s exit=%s\n" "$LINENO" "$rc" >&2; exit "$rc"' ERR
allocated=${SLURM_CPUS_PER_TASK:-}
if ! [[ "$allocated" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: submit through sbatch, not bash.\n' >&2
    exit 2
fi
# Reject an accidental portal override to a paid account/partition.
if [[ ${SLURM_JOB_ACCOUNT:-} != free || ${SLURM_JOB_PARTITION:-} != freecpuQ ]]; then
    printf 'ERROR: expected account=free partition=freecpuQ; got account=%s partition=%s.\n' \
        "${SLURM_JOB_ACCOUNT:-unset}" "${SLURM_JOB_PARTITION:-unset}" >&2
    exit 2
fi
project_dir=${SLURM_SUBMIT_DIR:?Submit with sbatch from the repository root}
cd "$project_dir"
clock_epoch() { date +%s.%N; }
job_start=$(clock_epoch)
shell_report="$project_dir/free-timing-shell-${SLURM_JOB_ID}.tsv"
printf 'Stage\tSeconds\n' > "$shell_report"
shell_stage() {
    local ended seconds
    ended=$(clock_epoch)
    seconds=$(awk -v end="$ended" -v start="$2" 'BEGIN {printf "%.6f",end-start}')
    printf '[SHELL TIME] %s %s s\n' "$1" "$seconds"
    printf '%s\t%s\n' "$1" "$seconds" >> "$shell_report"
}
printf 'ALLOCATION: job=%s account=%s partition=%s qos=%s cpus=%s mem_per_node_MB=%s\n' \
    "$SLURM_JOB_ID" "$SLURM_JOB_ACCOUNT" "$SLURM_JOB_PARTITION" \
    "${SLURM_JOB_QOS:-unknown}" "$allocated" "${SLURM_MEM_PER_NODE:-unknown}"
printf 'Workers are chosen by run_kemeny.m from enabled k values and actual CPU allocation.\n'
for required in run_kemeny.m src/kemeny_feedback.m third_party/bini/kemenydirect.m third_party/du/kemeny_k_parallel.m; do
    test -f "$required" || { printf 'ERROR: Missing %s\n' "$required" >&2; exit 2; }
done
if ! grep -q 'RUNNER_SCHEMA: multi-partitions-v2' run_kemeny.m; then
    printf 'ERROR: Upload the matching multi-partition run_kemeny.m together with this script.\n' >&2
    exit 2
fi
phase_start=$(clock_epoch)
module load MATLAB/R2021b
command -v matlab
shell_stage module_load "$phase_start"
phase_start=$(clock_epoch)
mkdir -p "$project_dir/runs"
run_dir="$project_dir/runs/free-${SLURM_JOB_ID}"
mkdir "$run_dir"
cp run_kemeny.m "$run_dir/"
cp -R src third_party matrix "$run_dir/"
export MATLAB_PREFDIR="$run_dir/matlab_preferences"
mkdir "$MATLAB_PREFDIR"
unset OMP_NUM_THREADS MKL_NUM_THREADS OPENBLAS_NUM_THREADS KEMENY_WORKERS
export KEMENY_OUTPUT_DIR="$run_dir"
cd "$run_dir"
shell_stage prepare_job_files "$phase_start"
printf 'RESULT_DIRECTORY: %s\n' "$run_dir"
cat > launch_matlab.sh <<'KEMENY_LAUNCH'
#!/bin/bash
set -eo pipefail
{
    hostname
    lscpu
    awk '/^Cpus_allowed_list:|^Mems_allowed_list:/' /proc/self/status
} > node-environment.txt
export KEMENY_MATLAB_LAUNCH_EPOCH=$(date +%s.%N)
set +e
/usr/bin/time -v -o process-resources.txt matlab -batch "run_kemeny('OutputDir',getenv('KEMENY_OUTPUT_DIR'))"
rc=$?
printf 'MATLAB exit code: %s\n' "$rc"
exit "$rc"
KEMENY_LAUNCH
phase_start=$(clock_epoch)
printf '[BEGIN SHELL] srun_and_MATLAB_TOTAL\n'
if srun --ntasks=1 --cpus-per-task="$allocated" bash ./launch_matlab.sh; then
    job_rc=0
else
    job_rc=$?
fi
shell_stage srun_and_MATLAB_TOTAL "$phase_start"
shell_stage job_shell_TOTAL "$job_start"
if [[ "$job_rc" -ne 0 ]]; then
    printf 'KEMENY_JOB_FAILED: exit=%s\n' "$job_rc" >&2
    exit "$job_rc"
fi
printf 'KEMENY_JOB_SUCCESS\n'
