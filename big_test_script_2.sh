#!/bin/bash

rm 1465/*/jobs/*/powers.*.rda
rm 1465/*/jobs/*/powers.one.*.bed
rm 1465/*/jobs/*/powers.two.*.bed

for f in 1465/power.1465-heart_BulkDNA_WGSb_job_scripts/*.sh; do
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA $f
done

while true
do
if [ $(bjobs | grep LIRA | wc -l) -eq 0 ]
then
break
fi
done

./lira collect_power -s config-example.txt -b config-example-bulk.txt
./lira bootstrap_germline -s config-example.txt -b config-example-bulk.txt
./lira call_ssnvs -s config-example.txt -b config-example-bulk.txt
