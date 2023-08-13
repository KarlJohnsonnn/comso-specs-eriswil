#!/bin/bash
#SBATCH --job-name=dummy_job
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --account=bb1262
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mem=0

# submit a dummy job
echo "Submitting a dummy job"
jobid=$(sbatch --parsable run_COSMO-SPECS_levante)
echo "Job submitted with ID: $jobid"

while true; do
    # Check the node that the job is assigned to
    assigned_node=$(squeue -j $jobid -o "%N" | tail -n +2)
    echo $assigned_node

    # Check if assigned_node matches the form "lXXXXX"
    if [[ $assigned_node =~ ^l[0-9]{5}$ ]]; then
        echo "Job is assigned to node: $assigned_node"
        break
    else
        echo "Job has not been assigned to a node yet, waiting... $assigned_node"
        sleep 5
    fi
done

# Print the job ID one more time at the end
echo "Job ID is: $jobid"


echo "cancel $jobid"
scancel $jobid
