# TXtools
Tool for processing 10X alignments


## Dependency

+ Python 2.7.10 or later

+ pysam 0.12 or later

## Usage

1. **Trim clipped reads**

	`python ./x_toolbox.py -T -b input.sorted.bam -r ${REF_FILE} -o output_trimmed.bam -p ./tmp -n 11`

2. **Convert reference indexed alignment to barcode indexed alignments**
	
	`python ./x_toolbox.py -C -b input.sorted.bam -o output_barcode_indexed.bam -p ./tmp -n 11`

3. **Convert barcode indexed alignments to reference based alignment**
	
	`python ./x_toolbox.py -K -b input.sorted.bam -i header_template_bam.bam -o output_barcode_indexed.bam`

4. **Haplotype block statistic**

	`python ./x_toolbox.py -B -b input.sorted.bam  -o block_info.txt -n 11 -p ./tmp/`

5. **Generate barocde matrix**

	`python ./x_toolbox.py -G -b input.sorted.bam -o out_matrix.bed -n 11 -p ./tmp/ -w 10000 -k 250`

4. **Merge alignments**
	
	`python ./x_toolbox.py -M -i merge_brain_mixed_list.txt -b ${BAM_FILE} -o out_merged.bam -n 11 -p ./tmp/`
