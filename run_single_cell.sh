#!/bin/bash

#SBATCH -n 10                               # 1 core
#SBATCH -t 2-15:00                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH --mem=30G
#SBATCH -p park                           # Run in short partition
#SBATCH -o hostname_%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=END                    # ALL email notification type
#SBATCH --mail-user=chong.simonchu@gmail.com  # Email to which notifications will be sent
#SBATCH --account=park_contrib

#./lira setup -c config-example-female.txt
#./lira split -c config-example-female.txt -p 
#./lira plink -c config-example-female.txt 
#./lira plink -c config-example-female.txt -a 5
#./lira compare -s config-example-female.txt -b config-example-bulk-female.txt -p -w 
#./lira ppower -s config-example-female.txt -b config-example-bulk-female.txt 
#./lira ppower -s config-example-female.txt -b config-example-bulk-female.txt -a 5
./lira varcall -s config-example-female.txt -b config-example-bulk-female.txt -o 

#./lira setup -c config-example-female.txt
#./lira split -c config-example-female.txt -p
#./lira plink  -c config-example-female.txt
#./lira compare -s config-example-female.txt -b config-example-bulk-female.txt -p -w -o
#./lira ppower  -s config-example-female.txt -b config-example-bulk-female.txt
#./lira varcall -s config-example-female.txt -b config-example-bulk-female.txt
#
#/n/data1/hms/dbmi/park/simon_chu/projects/LiRA2/lira_0514/4643_MDA_1/power.4643_Bulk-Liver_job_scripts/1_1.sh
#sh /n/data1/hms/dbmi/park/simon_chu/projects/LiRA2/lira_0514/4643_Bulk-Liver/job_scripts/2_15.sh
