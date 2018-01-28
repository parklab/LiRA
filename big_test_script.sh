#!/bin/bash
#rm -r 1465-2-C1 
#./lira setup -c config-example.txt
#for i in {1..22}
#do
#bsub -q short -W 12:0 -o bsub_out -e bsub_error 12:0 ./lira split_genome -c config-example.txt -m $i
#done
#./lira setup -c config-example.txt
#for i in {1..22}
#do
#bsub -q park_short -W 12:0 -o bsub_out -e bsub_error ./lira split_genome -c config-example.txt -m $i
#bsub -q park_short -W 12:0 -o bsub_out -e bsub_error ./lira split_genome -c config-example-bulk.txt -m $i
#bsub -q short -W 12:0 -o bsub_out -e bsub_error ./lira get_vcf_info -c config-example-bulk.txt -m $i
#bsub -q short -W 12:0 -o bsub_out -e bsub_error ./lira get_vcf_info -c config-example.txt -m $i
#done
#bsub -q short -W 12:0 -o bsub_out -e bsub_error ./lira split_genome -c config-example-bulk.txt -m 2

while read R; do
sed -i 's#^Rscript#checkR.sh;\nLD_LIBRARY_PATH=/tmp/R-3.3.1/lib /tmp/R-3.3.1/bin/Rscript#g' $R
echo "sleep 60" >> $R
bsub -q park_short -W 12:0 -o bsub_out -e bsub_error $R
done< <(cat tmp)

#for i in {1..22}
#do
#bsub -q park_short -W 12:0 -o bsub_out -e bsub_error ./lira get_alt_counts -c config-example-bulk.txt -m $i
#done
