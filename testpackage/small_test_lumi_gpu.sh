#!/bin/bash -l
#SBATCH --job-name=vlasi_g_tp
##SBATCH --partition=dev-g
#SBATCH --partition=small-g

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
##SBATCH --ntasks-per-node=8
##SBATCH --gpus-per-node=8

#SBATCH --time=24:00:00
##SBATCH --time=3:00:00
#SBATCH --account=project_462000358
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --hint=nomultithread

# If 1, the reference vlsv files are generated
# If 0 then we check the v1
create_verification_files=0

# Folder for all reference data 
reference_dir="/scratch/project_462000358/testpackage/"
cd $SLURM_SUBMIT_DIR

bin="/scratch/project_462000358/testpackage/vlasiator_gpu_2309_tp"
diffbin="/scratch/project_462000358/testpackage/vlsvdiff_DP_gpu"

# compare agains which revision?
reference_revision="current"

# place before exec
#LD_PRELOAD=/users/marbat/git/vlasiator-mempool/libpreload-me-2309.so

# set up GPU/CPU bindings
cat << EOF > select_gpu_${SLURM_JOB_ID}
#!/bin/bash
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
export OMP_NUM_THREADS=7
LD_PRELOAD=/users/marbat/git/vlasiator-mempool/libpreload-me-2309.so exec \$*
EOF
chmod +x ./select_gpu_${SLURM_JOB_ID}
# this should set the ordering correctly: "4 5 2 3 6 7 0 1"
CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

# module load LUMI/22.08
# module load partition/G
# module load cpeAMD
# module load rocm/5.3.3
# module load Boost/1.79.0-cpeAMD-22.08
module load LUMI/23.09
module load partition/G
module load cpeAMD
module load rocm/5.6.1
module load Boost/1.82.0-cpeAMD-23.09
module load papi/7.0.1.1
module load Eigen/3.4.0
module list

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_SUPPORT_ENABLED=1
# Turn off forced managed memory paging
export HSA_XNACK=0
# use extra threads for MPI in background
export MPICH_ASYNC_PROGRESS=1
# allow 16 in-parallel queues
export GPU_MAX_HW_QUEUES=16
# Command for running tests and diffs with MPI
# run_command="srun --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu "
# No MPI testing for now
run_command="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
small_run_command="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
run_command_tools="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID}"

umask 007

echo "Running on ${SLURM_NTASKS} mpi tasks with ${OMP_NUM_THREADS} threads per task on ${SLURM_JOB_NUM_NODES} nodes"

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
wait 2

# hello jobstep verification
# srun --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu ${SLURM_SUBMIT_DIR}/hello_jobstep

# cleanup
rm -f ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID}
