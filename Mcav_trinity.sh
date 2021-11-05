#!/bin/bash
#SBATCH -J mcav_trinity            # job name
#SBATCH -o mcav_trinity.o%j        # output and error file name (%j expands to jobID)
#SBATCH -e mcav_trinity.e%j        # name of stderr error file. 
#SBATCH -N 1                       # number of nodes requested
#SBATCH -n 32                      # total number of mpi tasks requested
#SBATCH -p nvdimm                  # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00                # run time (hh:mm:ss) - 48 hours
#SBATCH -A Immunity-to-Communit    #charge job to the specified project 
#SBATCH --mail-user=kelsey.beavers@uta.edu
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

# Other commands must follow all #SBATCH directives...
module load tacc-singularity

#set the path
PATH=$PATH:/scratch/06825/tg861249

# run the executable
singularity exec -e trinityrnaseq.v2.12.0.simg Trinity --normalize_reads --seqType fq --grid_node_CPU 21 --grid_node_max_memory 300G --max_memory 300G --SS_lib_type FR --left UVI_11_trim_1.fq.gz,UVI_13_trim_1.fq.gz,UVI_15_trim_1.fq.gz,UVI_17_trim_1.fq.gz,UVI_19_trim_1.fq.gz --right UVI_11_trim_2.fq.gz,UVI_13_trim_2.fq.gz,UVI_15_trim_2.fq.gz,UVI_17_trim_2.fq.gz,UVI_19_trim_2.fq.gz --CPU 25