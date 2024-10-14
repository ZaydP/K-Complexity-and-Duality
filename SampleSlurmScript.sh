#!/bin/sh
#SBATCH --account=maths
#SBATCH --partition=ada
#BATCH --time=02:30:00
#SBATCH --nodes=1 --ntasks=30
#SBATCH --job-name="K-Complexity_test_script"
#SBATCH --mail-user=pndzay001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL

module load software/mathematica-12.3

wolframscript -file K-Complexity_test_script.wls
