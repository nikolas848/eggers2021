#! /bin/bash

# findPeaks.sbatch

#SBATCH -J meme   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 1                         # 1 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-2:00                    # Running time of 2 hours
#SBATCH --mem 8000                  # Memory request
#SBATCH -o out/meme_%A.out          # Standard output
#SBATCH -e out/meme_%A.err          # Standard error

# grab out filename

module load meme/5.0.2 

meme ${f} -oc ${f}memeout -mod zoops -dna -revcomp -nmotifs 1

meme ${f} -oc ${f}memeoutanr -mod anr -dna -revcomp -nmotifs 1