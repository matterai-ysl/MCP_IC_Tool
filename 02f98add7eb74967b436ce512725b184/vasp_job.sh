#!/bin/bash
#SBATCH --job-name=02f98add7eb74967b436ce512725b184
#SBATCH --partition=p1
#SBATCH -N 2
#SBATCH -n 56
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load vasp/6.3.2-intel
source /data/app/intel/oneapi-2023.2/setvars.sh >/dev/null 2>&1
ulimit -s unlimited
ulimit -l unlimited

echo "=== 作业信息 ==="
echo "作业ID: $SLURM_JOB_ID"
echo "分区: $SLURM_JOB_PARTITION"
echo "节点数: $SLURM_JOB_NUM_NODES"
echo "总任务数: $SLURM_NPROCS"
echo "每节点任务数: $SLURM_NTASKS_PER_NODE"
echo "节点列表: $SLURM_JOB_NODELIST"

echo "=== 开始VASP计算 ==="
mpirun -np $SLURM_NPROCS vasp_std>result.log 2>&1
echo "VASP计算完成
        