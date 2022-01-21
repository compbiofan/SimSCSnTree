#!/bin/bash
#SBATCH --job-name=run.0
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=xf2@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/xf2/github/SingleCellCNABenchmark/evaluation/run.0.out
#SBATCH --error=/home/xf2/github/SingleCellCNABenchmark/evaluation/run.0.err
#python count_overlap.py -f /projects/nakhleh/GB20/p1/rep1/newoverlappingCNA.csv -i /projects/nakhleh/GB20/p1/rep1/leaves.txt -t /projects/nakhleh/GB20/p1/rep1/from_first_step.tree.npy
python count_overlap.py -f /projects/nakhleh/GB20/p1/rep1/newoverlappingCNA.csv -i /projects/nakhleh/GB20/p1/rep1/leaves.txt -t /projects/nakhleh/GB20/p1/rep1/from_first_step.tree.npy -a /projects/nakhleh/xf2/reference/hg19.fa.fai -T 0.33
