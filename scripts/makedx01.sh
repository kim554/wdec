#!/bin/bash
# submit_sequential_models.sh (in each folder, e.g., folder01/submit_sequential_models.sh)

# Configuration for this specific folder
EXECUTABLE_NAME="makedx01" # <<< IMPORTANT: CHANGE THIS FOR EACH FOLDER (e.g., makedx02, makedx03, etc.)
OUTPUT_DIR="slurm_job_files"
main_param_file="gridparameters.torun" # The file containing your parameters

# --- SLURM Configuration for this single job ---
# Set these as appropriate for your cluster and job characteristics
SLURM_PARTITION="nodes"    # Example: "debug", "batch", "nodes" - REQUIRED for most clusters
SLURM_TIME="7-00:00:00"    # D-HH:MM:SS - Set a sufficiently long time limit for ALL 171 models
                           # 171 models * 30 seconds/model = 5130 seconds = ~1.4 hours
                           # So, "0-02:00:00" (2 hours) or "0-03:00:00" (3 hours) for safety is reasonable.
SLURM_MEM="4G"            # Memory for the entire job (if makedx is lightweight, 16G is probably ample for one instance)
SLURM_CPUS_PER_TASK=1      # Since makedx needs to run sequentially per folder, 1 CPU is sufficient
# --- End SLURM Configuration ---

# Ensure the output directory exists
mkdir -p "${OUTPUT_DIR}"

# Add this line to clear its contents at the start of each run
# This ensures you're starting with a clean slate for the test
rm -f "${OUTPUT_DIR}"/*

# Calculate the total number of parameter lines
np=$(wc -l < "$main_param_file")

if [ "$np" -eq 0 ]; then
    echo "Error: The file '${main_param_file}' is empty or does not exist." >&2
    exit 1
fi

echo "Found ${np} parameter lines in '${main_param_file}'. Submitting a single job to run them sequentially."

# Construct the SLURM batch script for this single, long-running job
SLURM_JOB_SCRIPT_NAME="sequential_job_${EXECUTABLE_NAME}.slurm"
SLURM_JOB_SCRIPT="${OUTPUT_DIR}/${SLURM_JOB_SCRIPT_NAME}"

# ... (rest of your submit_sequential_models.sh script before cat << EOF) ...

cat << EOF > "${SLURM_JOB_SCRIPT}"
#!/bin/bash
#
#SBATCH -n 1
#SBATCH --cpus-per-task=${SLURM_CPUS_PER_TASK}
#SBATCH --mem=${SLURM_MEM}
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -t ${SLURM_TIME}
#SBATCH -J ${EXECUTABLE_NAME}
#SBATCH -o ${OUTPUT_DIR}/slurm_main_output.log # OUTPUT_DIR expanded by outer script
#SBATCH -e ${OUTPUT_DIR}/slurm_main_error.log # OUTPUT_DIR expanded by outer script

# Export necessary environment variables
export MESA_DIR="/home/axk55/mesa-r22.11.1"
export EXECUTABLE_PATH="\$(pwd)/${EXECUTABLE_NAME}" # \$ for expansion by SLURM job
export MAIN_PARAM_FILE="\$(pwd)/${main_param_file}" # \$ for expansion by SLURM job

echo "--- ulimit settings for job \${SLURM_JOB_NAME} (Job ID \${SLURM_JOB_ID}) ---"
ulimit -a
echo "---------------------------------------------------"

echo "Job_ID \${SLURM_JOB_ID} starting on \$(hostname)"
echo "Using MESA_DIR: \${MESA_DIR}"
echo "Total parameter lines to process: ${np}"

# Ensure we are in the original submission directory where this job was launched
cd "\${SLURM_SUBMIT_DIR}" || { echo "Error: Could not change to SLURM_SUBMIT_DIR \${SLURM_SUBMIT_DIR}" >&2; exit 1; }

# --- Diagnostic: Enable shell debugging (VERY VERBOSE! Temporary for debugging) ---
set -x

# --- THE SEQUENTIAL LOOP ---
for i in \$(seq 1 ${np}); do # \$ for expansion by SLURM job
    JOB_ITERATION_NUM=\$i # Use a new variable for clarity inside the loop, \$ for expansion by SLURM job

    echo "--- Processing model \$i / ${np} ---"
    echo "Current working directory: \$(pwd)"

    # Extract the specific parameter line for *this* model
    line_content_for_this_model=\$(sed -n "\${JOB_ITERATION_NUM}p" "\${MAIN_PARAM_FILE}")

    if [ -z "\${line_content_for_this_model}" ]; then
        echo "Error: Parameter line \${JOB_ITERATION_NUM} from \${MAIN_PARAM_FILE} is empty or missing. Skipping this model." >&2
        continue
    fi

    echo "\${line_content_for_this_model}" > gridparameters
    echo "Files in current directory BEFORE makedx run:"
    ls -lha .

    echo "Executing: \"\${EXECUTABLE_PATH}\""
    # Fixed stdout/stderr logs to use JOB_ITERATION_NUM
    "\${EXECUTABLE_PATH}" > "${OUTPUT_DIR}/makedx_stdout_model\${JOB_ITERATION_NUM}.log" 2> "${OUTPUT_DIR}/makedx_stderr_model\${JOB_ITERATION_NUM}.log"

    exit_status=\$?
    echo "makedx executable finished with exit status: \${exit_status}"
    echo "Files in current directory AFTER makedx run:"
    ls -lha .
    
    if [ \${exit_status} -eq 0 ]; then
        echo "Execution successful for \${SLURM_JOB_NAME} (Model \${JOB_ITERATION_NUM})."
        
        # Check for expected output files and move them
        if [ -f calcperiods ]; then
            echo "Found calcperiods. Moving to ${OUTPUT_DIR}/model\${JOB_ITERATION_NUM}.calcperiods"
            # --- THE ACTUAL FIX FOR `mv` DESTINATION ---
            mv calcperiods "${OUTPUT_DIR}/model\${JOB_ITERATION_NUM}.calcperiods"
            mv_calcperiods_status=\$?
            echo "mv calcperiods status: \${mv_calcperiods_status}"
        else
            echo "WARNING: calcperiods not found for model \${JOB_ITERATION_NUM}!" >&2
        fi

        if [ -f radii.dat ]; then
            echo "Found radii.dat. Moving to ${OUTPUT_DIR}/model\${JOB_ITERATION_NUM}.radii.dat"
            # --- THE ACTUAL FIX FOR `mv` DESTINATION ---
            mv radii.dat "${OUTPUT_DIR}/model\${JOB_ITERATION_NUM}.radii.dat"
            mv_radii_status=\$?
            echo "mv radii.dat status: \${mv_radii_status}"
        else
            echo "WARNING: radii.dat not found for model \${JOB_ITERATION_NUM}!" >&2
        fi

        rm -f gridparameters
    else
        echo "Error: \$(basename "\${EXECUTABLE_PATH}") failed for \${SLURM_JOB_NAME} (Model \${JOB_ITERATION_NUM}) with exit code \${exit_status}" >&2
    fi
    echo "Files in current directory AFTER cleanup and move:"
    ls -lha .
    echo "--- End processing model \$i ---"
done
# --- END SEQUENTIAL LOOP ---

echo "Job_ID \${SLURM_JOB_ID} finished all sequential models."

# --- Diagnostic: Disable shell debugging ---
set +x

# Exit with success if all models processed, otherwise with error
if [ \$(grep -c "Error: .* failed" "${OUTPUT_DIR}/slurm_main_error.log") -eq 0 ]; then
    exit 0
else
    exit 1
fi
EOF

# Submit the single job
echo "Submitting single job script: ${SLURM_JOB_SCRIPT}..."
job_id=$(sbatch --parsable "${SLURM_JOB_SCRIPT}")
# Clean up any whitespace (including newlines) from the job ID
job_id=$(echo "${job_id}" | tr -d '[:space:]') 

# --- Add this line temporarily for debugging if it still fails ---
# echo "DEBUG: Cleaned job_id content: ---${job_id}---"

# Check if sbatch command was successful (this check uses $? from the sbatch command itself)
if [ $? -ne 0 ]; then # This $? refers to the exit status of 'sbatch'
    echo "Error submitting job. Slurm command failed." >&2
    cat "${SLURM_JOB_SCRIPT}" >&2 # Output the script content for debugging
    exit 1
fi

# This check ensures the captured job_id is purely numeric after cleaning
if [[ "$job_id" =~ ^[0-9]+$ ]]; then
    echo "Successfully submitted single job (ID: ${job_id}) to run ${np} models sequentially."
    echo "Monitor with 'squeue -u \$USER' or 'squeue -j ${job_id}'"
    echo "Main output will be in ${OUTPUT_DIR}/slurm_main_output.log"
    echo "Individual model outputs will be in ${OUTPUT_DIR}/makedx_stdout_modelX.log, etc."
else
    # This should now ideally only trigger if sbatch somehow failed to return a number at all,
    # or if the cleaning failed for an unexpected reason.
    echo "Error: sbatch did not return a valid job ID after cleaning. Output: '${job_id}'" >&2
    echo "Please check Slurm configuration and permissions." >&2
    exit 1
fi

echo "All setup for this folder complete."
