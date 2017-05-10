#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv)
{
  int rank;
  char hostname[256];
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  gethostname(hostname,255);
  printf("Hello world!  I am process number: %d on host %s\n", rank, hostname);
  MPI_Finalize();
  return 0;
}
// module load compiler/gcc/4.9.3/compilervars
//module load compiler/mpi/mpich/3.2/gnu
#PBS -l select=10:mem=10GB:ncpus=12:ngpus=1:K20GPU=false
### Specify "wallclock time" required for this job, hhh:mm:ss
### Specify output error files
#PBS -o logs.s.out
#PBS -e logs.s.err
#PBS -l walltime=00:05:00
#### Get environment variables from submitting shell
#PBS -V
#PBS -l software=CUDA
# After job starts, must goto working directory.
# $PBS_O_WORKDIR is the directory from where the job is fired.
echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo $PBS_NODEFILE
echo "==============================="
cd $PBS_O_WORKDIR
module load compiler/gcc/4.9.3/compilervars
module load mpi/openmpi/1.10.0/gcc/mpivars
#job
mpirun -n 10 --hostfile $PBS_NODEFILE --mca mpi_preconnect_mpi 1 /home/cse/btech/cs1130203/pmlib/testSuite/matrixMultiply/build/linux/release/matrixMultiply.exe > /home/cse/btech/cs1130203/log_uni.txt
