#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N hoagie_dynamix-ga
#$ -m bae -M andyras@gmail.com
#$ -o stdout.log
#$ -e stderr.log
#$ -pe orte 128

function onKill {
echo "DANGER!"
echo "DANGER!"
echo "Job terminated unexpectedly!"
echo "JOB END TIME: $(date +"%F %T")"
echo "JOB DURATION: $(($(date +%s) - ${startTime})) s"
cd ${OWD}
exit
}

trap onKill 2 9 15

# Print debug info about job
echo "# Job information"
echo "HOSTNAME: $(hostname)"
echo "DIRECTORY: $(pwd)"
echo "JOB START TIME: $(date +"%F %T")"
startTime=$(date +%s)
echo ""

# create scratch directory
echo "create scratch directory"
SCRATCH_JOB_DIR=${TMPDIR}/$(basename ${SGE_CWD_PATH})
cp -rf ${SGE_CWD_PATH} ${SCRATCH_JOB_DIR}
rm -f ${SCRATCH_JOB_DIR}/std{err,out}.log

# Go to SGE working directory
echo "Going to SGE working directory: ${SCRATCH_JOB_DIR}"
cd ${SCRATCH_JOB_DIR}

# Run program
echo "Running program"
echo ""

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

. $HOME/scripts/module_GAlib-mpi load
OPENMPI=/opt/openmpi/1.6.5/intel/14.0.2
${OPENMPI}/bin/mpirun --mca ras_gridengine_show_jobid 1 --mca ras_gridengine_debug 1 --mca ras_gridengine_verbose 30 --mca ras_base_verbose 30 --mca btl_base_verbose 30 --display-allocation --prefix ${OPENMPI} -x LD_LIBRARY_PATH ${SGE_CWD_PATH}/bin/dynamix-ga
#mpirun ${SGE_CWD_PATH}/bin/dynamix-ga

echo ""
echo "JOB END TIME: $(date +"%F %T")"
echo "JOB DURATION: $(($(date +%s) - ${startTime})) s"
cd ${OWD}

# copy job dir contents back to where they started
echo "copying from scratch to original directory"
cp -rf ${SCRATCH_JOB_DIR}/* ${SGE_CWD_PATH}/
cd $SGE_CWD_PATH
