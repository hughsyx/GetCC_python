#!/bin/bash
#
################
#set a job name
#SBATCH --job-name=CC_mlt2
################
#a file for job output
#SBATCH --output=cc_mlt_15.out2
################
#a file for errors from the job
#SBATCH --error=cc_mlt_15.err
################
#time I need hh:mm:ss
#SBATCH --time=48:00:00
################
#quality of service
#SBATCH --qos=normal
################
#number of nodes I need
#SBATCH --nodes=1
################
#memory per node
#SBATCH --mem=4000
################
#tasks to run per node
#SBATCH --ntasks-per-node=16
################

python /scratch/users/yixiao2/codes/python/run_GetANIR_para_mult_rcv.py YN.BCCC
