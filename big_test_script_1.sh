#!/bin/bash
rm -r 1465-2-C1
rm -r 1465-heart_BulkDNA_WGSb 
./lira setup -c config-example.txt
./lira setup -c config-example-bulk.txt

for i in {1..22}
do
bsub -q park_short -J LIRA -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira split_genome -c config-example.txt -m $i
bsub -q park_short -J LIRA -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira split_genome -c config-example-bulk.txt -m $i
bsub -q park_short -J LIRA -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira get_vcf_info -c config-example-bulk.txt -m $i
bsub -q park_short -J LIRA -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira get_vcf_info -c config-example.txt -m $i
done

while true
do
if [ $(bjobs | grep LIRA | wc -l) -eq 0 ]
break
fi
done

while read R; do
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA 1465-2-C1/job_scripts/$R
done< <(ls 1465-2-C1/job_scripts)

while read R; do
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA 1465-heart_BulkDNA_WGSb/job_scripts/$R
done< <(ls 1465-heart_BulkDNA_WGSb/job_scripts)

for i in {1..22}
do
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira get_alt_counts -c config-example-bulk.txt -m $i
done

while true
do
if [ $(bjobs | grep LIRA | wc -l) -eq 0 ]
break
fi
done

for i in {1..22}
do
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira collect_results -c config-example.txt -m $i
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira collect_results -c config-example-bulk.txt -m $i
done

while true
do
if [ $(bjobs | grep LIRA | wc -l) -eq 0 ]
break
fi
done

for i in {1..22}
do
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error -J LIRA ./lira compare_to_bulk -s config-example.txt -b config-example-bulk.txt -m $i
done
