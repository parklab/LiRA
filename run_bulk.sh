#!/bin/bash

#SBATCH -n 12                               # 1 core
#SBATCH -t 2-15:00                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH --mem=60G
#SBATCH -p park                           # Run in short partition
#SBATCH -o hostname_%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=END                    # ALL email notification type
#SBATCH --mail-user=chong.simonchu@gmail.com  # Email to which notifications will be sent
#SBATCH --account=park_contrib


./lira setup -c config-example-bulk-female.txt
./lira split -c config-example-bulk-female.txt -p
./lira plink -c config-example-bulk-female.txt
./lira plink -a 5 -c config-example-bulk-female.txt

#./lira setup -c config-example-female.txt
#./lira split -c config-example-female.txt -p
#./lira plink -c config-example-female.txt
#./lira compare -s config-example-female.txt -b config-example-bulk-female.txt -p -w
#./lira ppower -s config-example-female.txt -b config-example-bulk-female.txt
#./lira varcall -s config-example-female.txt -b config-example-bulk-female.txt
