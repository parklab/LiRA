#!/bin/sh

#BSUB -n 1                  # Job will use 6 cores
#BSUB -W 100:00                 # Job will be killed if it runs longer than 3 hours
#BSUB -J LiRA_1465       # Job name
#BSUB -o %J.out             # Output file name. %J is automatically replaced by the job ID.
#BSUB -e %J.err             # Error file name. %J is automatically replaced by the job ID.
#BSUB -q park_medium             # Queue name.
#BSUB -R "rusage[mem=5000]" # Use 1GB memory PER CORE, or 8GB total

sh big_test_script_1.sh
