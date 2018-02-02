#!/bin/bash

#SBATCH -n 1                               # 1 core
#SBATCH -t 5-00:00                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH --mem=100G
#SBATCH -p priority                           # Run in short partition
#SBATCH -o hostname_%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=END                    # ALL email notification type
#SBATCH --mail-user=chong.simonchu@gmail.com  # Email to which notifications will be sent

#./lira split_genome  --config config-example.txt --chr 1
#./lira get_vcf_info --config config-example.txt --chr 1
./lira get_alt_counts --config config-example-bulk.txt --chr 1
/n/data1/hms/dbmi/park/simon_chu/projects/LiRA/test/1465/job_scripts/1_1.sh
./lira check_results --config config-example.txt
./lira collect_results --config config-example.txt --chr 1
./lira compare_to_bulk --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --chr 1  
./lira collect_compare --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
./lira power_jobs --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
