# LiRA
LiRA ("Linked Read Analysis") is a computational tool for removing amplification artifacts from single-cell DNA sequencing data and estimating mutation rates in single cells. LiRA utilizes the fact that true somatic variation appears in linkage with known, germline heterozygous variation, while artifacts do not.  It uses read-backed phasing to evaluate candidate somatic single-nucleotide variants (sSNVs) occurring near population-polymorphic germline heterozygous sites (gHets), and from the number of such regions covered across the genome infers sensitivity and estimates the genome-wide rate of somatic mutation.  

LiRA's sensitivity is limited to regions of the genome close enough to gHets to permit read-backed phasing, but achieves unparalleled specificity.  In contrast to previous methods for single-cell analysis, LiRA can operate robustly on variants detected in only one cell, and as such can be used to analyze mutations and mutation rates in arbitrarily small cell populations and post-mitotic cells.

For more information, please see our [bioarchive paper](https://doi.org/10.1101/211169) and [our related paper in science](https://doi.org/10.1126/science.aao4426) examining variation in mutational burden in single neurons across age and tissue type, and in progeroid disorders.

## Requirements
The following include tested versions in parenthensis; later versions are likely to still work.
0. linux / macOS
1. python (2.7.12)
	+ argparse
	+ os
2. samtools (1.2)
3. bcftools (1.2)
4. htslib (1.2.1)
5. bedtools (2.23.0)
6. picard (2.7.1)
6. R (3.3.1)
	+ stringr
	+ digest
7. SHAPEIT2
8. SnpSift (optional, for annotating common variants)

Download [samtools, bcftools, and htslib](https://github.com/samtools)

Download [bedtools](https://github.com/arq5x/bedtools2/releases)

Download [SHAPEIT2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download)

Download [SnpSift](http://snpeff.sourceforge.net/download.html)

Download [picard](https://github.com/broadinstitute/picard/releases/tag/2.17.8)

## Installation
1. Clone this repository `git clone https://github.com/parklab/LiRA`
2. Set the `LIRA_DIR` environmental variable to the absolute path to your local copy of this repository (e.g. if it is in /home/me/LiRA, then add `export LIRA_DIR=/home/me/LiRA` to your `.bash_profile`)
3. Ensure that samtools, bcftools, and bedtools are in your `PATH`.
4. Set the `PICARD` environmental variable to the absolute path to your picard jar file.

## Inputs to LiRA (required)
  * A set of called variants in vcf format.  This callset should include both germline and prospective somatic variants, and must include at least one single-cell and bulk sample.  At the moment, to run LiRA the input vcf must have been generated using GATK haplotype caller, and must be annotated with 1000 genomes common allele frequencies in the info field 'DBSNP_CAF.'  In future versions, a LiRA tool will assist with this.  For now, this can be accomplished using the annotate tool from SnpSift and a dbSNP release from the NCBI database.  See: https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/, e.g. `common_all_20170710.vcf.gz`
  * The BAM files used to produce the vcf (one for the single cell, and one for bulk).
  * A phased vcf file with information on all common variant sites called in the bulk sample, obtained using SHAPEIT2. A common variant is one that has a defined DBSNP_CAF value (in other words, an annotated population allele frequency).

## Usage
You can set up a new analysis with LiRA by creating and editting new config files (see examples).  `<analysis_path>` and `<name>` (below) are both assigned values in their associated config files. LiRA commands amount to executing successive steps in the LiRA bioinformatic workflow, which must proceed in this order:

Note: as the examples provided here require large files not a part of this repository, they will not run; they are merely meant to illustrate syntax.

Steps 3a, 3b, and 3c can be run concurrently.

1. setup (single cell and bulk)
  * Arguments: Config file
  * Creates a new analysis directory and soft links to input files.
  * Example: 
```
./lira setup --config config-example.txt
```

2. split_genome (single cell and bulk)
  * Arguments: Config file, chromosome
  * Analyzes the genome-wide coverage profile of the input BAM file and chunks the genome into smaller pieces for analysis.  Bash scripts are created in the ANALYSIS_PATH/job_scripts directory. Must be run for all diploid chromosomes (1-22 for males, 1-22 + X for females).
  * Example: 
```
		./lira split_genome  --config config-example.txt --chr 1
		./lira split_genome  --config config-example.txt --chr 2
		./lira split_genome  --config config-example.txt --chr 3
```

3a. get_vcf_info (single cell and bulk)
  * Arguments: Config file, chromosome
  * Extracts relevant info from the input vcf file over the input chromosome.  Must be run for all diploid chromosomes.
  * Example: 
```
		./lira get_vcf_info --config config-example.txt --chr 1
		./lira get_vcf_info --config config-example.txt --chr 2
		./lira get_vcf_info --config config-example.txt --chr 3
```

3b. get_alt_counts (only bulk)
  * Arguments: Config file, chromosome
  * Gets raw counts of alternate allele coverage at all variant sites in the input vcf file. Must be run for all diploid chromosomes.
get_alt_counts (only necessary over bulk samples)
  * Example: 
```
		./lira get_alt_counts --config config-example-bulk.txt --chr 1
		./lira get_alt_counts --config config-example-bulk.txt --chr 2
		./lira get_alt_counts --config config-example-bulk.txt --chr 3
```

3c. (Outside LiRA) Run scripts produced by split_genome (single cell and bulk)
  * Each script analyzes linkage between variants over a subset of the genome.
  * Example:
```
		<analysis_path>/job_scripts/1_1.sh
		<analysis_path>/job_scripts/1_2.sh
		<analysis_path>/job_scripts/1_3.sh
```
  * There are a large number of these scripts, and it is intended that you use LiRA's built-in 'chunking' in conjunction with your own cluster job scheduler.

4. check_results (single cell and bulk)
  * Arguments: Config file
  * Verifies completion of the scripts from 3c; prints incomplete jobs and exits or verifies completeness.
  * Example:
```
./lira check_results --config config-example.txt
```

5. collect_results (single cell and bulk)
  * Arguments: Config file
  * Collects results from 3c into files for each diploid chromosome. Must be run for all diploid chromosomes.
  * Example: 
```
		./lira collect_results --config config-example.txt --chr 1
		./lira collect_results --config config-example.txt --chr 2
		./lira collect_results --config config-example.txt --chr 3
```

6. compare_to_bulk
  * Arguments: Single cell config file, bulk config file, chromosome
  * Compares linked variants in single cell and bulk data over a single diploid chromosome. Must be run for all diploid chromosomes.
  * Example: 
```
		./lira compare_to_bulk --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --chr 1      
		./lira compare_to_bulk --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --chr 2
		./lira compare_to_bulk --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --chr 3
```

7. collect_compare 
  * Arguments: Single cell config file, bulk config file
  * Collects the results of (6)
  * Example:
```
./lira collect_compare --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
```

8. power_jobs
  * Arguments: Single cell config file, bulk config file
  * Creates bash scripts to measure power in a chunked fashion across the genome.  Scripts are placed in <analysis_path (single-cell)>/power.<name (bulk)>_job_scripts
  * Example: 
```
./lira power_jobs --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
```

9. (Outside LiRA) Run scripts produced by power_jobs
  * Each script analyzes the power LiRA has to detect somatic variation over different regions of the genome.
  * Example:
```
		<analysis_path (single-cell)>/power.<name (bulk)>_job_scripts/1_1.sh
		<analysis_path (single-cell)>/power.<name (bulk)>_job_scripts/1_2.sh
		<analysis_path (single-cell)>/power.<name (bulk)>_job_scripts/1_3.sh
```
  * There are a large number of these scripts, and it is intended that you use LiRA's built-in 'chunking' in conjunction with your own cluster job scheduler.

10. check_power
  * Arguments: Single cell config file, bulk config file
  * Checks the completeness of (9); prints incomplete jobs and exits or verifies completeness.
  * Example: 
```
./lira check_power --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
```

11. collect_power
  * Arguments: Single cell config file, bulk config file
  * Collect the results of (9)
  * Example: 
```
./lira collect_power --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
```

12. bootstrap_germline
  * Arguments: Single cell config file, bulk config file
  * Select random, linked distance and size-matched germline variant sets to calibrate the relationshipbetween false positive rate and composite coverage.
  * Example: 
```
./lira bootstrap_germline --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
```

13. call_ssnvs
  * Arguments: Single cell config file, bulk config file
  * Call sSNVs and calculate genome-wide mutation rate using the results of compare_to_bulk and bootstrap_germline.  Produces an rda file (ssnvs.<name (bulk)>.rda) containing info about all somatic SNPs.  Also produces a rate plot (rate-plot.<name (bulk)>.rda) showing the sSNV rate vs. composite. coverage and the threshold decision.
  * Example: 
```
./lira call_ssnvs --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
```
