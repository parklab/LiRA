# LiRA
LiRA (**Li**nked **R**ead **A**nalysis) is a computational tool for removing amplification artifacts from single-cell DNA sequencing data and estimating mutation rates in single cells.  It uses read-backed phasing to evaluate candidate somatic single-nucleotide variants (sSNVs) occurring near population-polymorphic germline heterozygous sites (gHets), and from the number of such regions covered across the genome infers sensitivity and estimates the genome-wide rate of somatic mutation.  To detect and remove artifacts, LiRA utilizes the fact that true somatic variation appears in linkage with known, germline heterozygous variation, while artifacts do not.

While analysis using LiRA is limited to regions of the genome close to gHets, its specificity is unparalleled.  In contrast to previous methods for single-cell analysis, LiRA can operate robustly on variants detected in only one cell, and as such can be used to analyze mutations and mutation rates in arbitrarily small cell populations and post-mitotic cells.

For more information, please see our [paper in nature genetics](https://doi.org/10.1038/s41588-019-0366-2) and [our related paper in science](https://doi.org/10.1126/science.aao4426) examining variation in mutational burden in single neurons across age and tissue type, and in progeroid disorders.

## Requirements
The following include tested versions in parenthesis when applicable; later versions are likely to still work.  These instructions are designed to enable use of LiRA on human sequencing data aligned to NCBI build 37 (hg19).  While it is possible to analyze other organisms or other genome builds, doing so will require an appropriate haplotype set alternative and set of common SNPs (as opposed to the human polymorphisms provided in the DBSNP release and the haplotype sets provided by the 1000 genomes phase 3 release, respectively).

0. linux / macOS
1. python (2.7.12)
	+ argparse
	+ os
2. samtools (1.2)
3. bcftools (1.2)
4. htslib (1.2.1)
5. bedtools (2.23.0)
6. java (1.8)
7. picard (2.7.1)
8. R (3.3.1)
	+ stringr
	+ digest
9. SHAPEIT2 (Note: Ensure that downloaded `.txt`, `.hap.gz`, and `.legend.gz` files for chromosome X are in the same directory as the other downloaded files.) or EAGLE
10. SnpSift
11. A DBSNP VCF release
12. The 1000 genomes phase 3 integrated haplotype set
13. (Highly preferable) Access to cluster computing

Download [samtools, bcftools, and htslib](https://github.com/samtools)

Download [bedtools](https://github.com/arq5x/bedtools2/releases)

Download [SHAPEIT2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download) or [EAGLE](https://data.broadinstitute.org/alkesgroup/Eagle/)

Download [SnpSift](http://snpeff.sourceforge.net/download.html)

Download [picard](https://github.com/broadinstitute/picard/releases/tag/2.17.8)

Download [a DBSNP release](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/), e.g. `common_all_20170710.vcf.gz`

Download [the 1000 genomes phase 3 haplotype set](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html)

## Installation
1. Clone this repository `git clone https://github.com/parklab/LiRA`
2. Set the `LIRA_DIR` environmental variable to the absolute path to your local copy of this repository (e.g. if it is in `/home/me/LiRA`, then add `export LIRA_DIR=/home/me/LiRA` to your `.bash_profile`)
3. Ensure that `samtools`, `bcftools`, `bedtools`, `tabix`, `java`, `R`, and `shapeit` (if applicable) are in your `PATH`.
4. Edit the `global_config.txt` file to point the `SNPEFF`, `DBSNP`, `KGEN`, `PICARD`, and (if applicable) `EAGLE`, `EAGLE_HG19_REF`, and `EAGLE_HG38_REF` variables to the appropriate locations.
5. (Optional) Ensure that `lira` is in your `PATH`.

## Inputs to LiRA (required)
  * A set of called variants in VCF format.  This call set should include both germline and prospective somatic variants and must include at least one single-cell and bulk sample (but can include more).In other words, the single-cell and bulk samples must come from the VCF file. Running LiRA on VCF files generated over more samples will increase runtime but will provide a (potential) advantage, as it will enable viewing in all samples the evidence for all variation discovered across all samples.  At the moment, to run LiRA the input VCF must have been generated using GATK Haplotype Caller.  The variant quality score recalibration (VQSR) step in the GATK best practices workflow is unnecessary.
  * The BAM files used to produce the VCF (one for single cell, and one for bulk).

## Outputs of LiRA
  * A VCF file containing information about all germline and candidate somatic sites considered by LiRA.  This file is created in the 'varcall' step outlined below, and is output at the location `ANALYSIS_PATH/out.NAME[bulk].vcf.gz`.
  * A text file describing the inferred genome-wide rate of somatic mutation, output at the location `ANALYSIS_PATH/summary.NAME(bulk).txt`.
  * Bed files describing where LiRA had power to call somatic mutations.  LiRA outputs two bed files describing power in an allele-specific manner (one for each haploid genome copy at the location `ANALYSIS_PATH/powers.(one|two).NAME[bulk].bed`), and one bed file `ANALYSIS_PATH/powers.NAME[bulk].bed` describing the union of these two files and, in the fourth column, whether LiRA had power over one (1) or both (2) alleles.  
  * Log files showing what various commands are doing (various subdrectoreies under `ANALYSIS_PATH`, always called `log.txt`)

## Parallelization (recommended)
LiRA is computationally intensive, but 40X WGS data should finish in 1-2 days given access to cluster computing resources.  Provided with LiRA in `scripts/slurm.R` are two functions: `check.jobs` and `submit.jobs`.  These scripts manage parallelization using slurm.  If your cluster resources have a different job management system, running LiRA in parallel is still possible: define new `check.jobs` and `submit.jobs` functions in an Rscript and LiRA will be able to use them.  The required behavior of these functions is described in `scripts/slurm.R`. After writing your own script, move it to `$LIRA_DIR/scripts` and point LiRA to it by editing the `PARALLEL_SCRIPT` parameter in your local copy of `global-config.txt`.  If you must run LiRA serially, you can alter `PARALLEL_SCRIPT` to point to the provided Rscript `serial.R`.

## Usage
You can set up a new analysis with LiRA by creating and editing new config files (see examples).  `<analysis_path>` and `<name>` (below) are both assigned values in their associated config files. LiRA commands amount to executing successive steps in the LiRA bioinformatic workflow, which must proceed in this order:

Note: as the examples provided here require large files not a part of this repository, they will not run; they are merely meant to illustrate syntax.  Also, analysis of 10x genomics synthetic long read data is still being implemented and as such is currently unsupported.

1. setup (single cell and bulk)
  * Arguments: Config file
  * Creates a new analysis directory and soft links to input files.
  * Example: 
```
    ./lira setup --config config-example.txt
```

2. split (single cell and bulk)
  * Arguments: Config file, chromosome
  * Analyzes the genome-wide coverage profile of the input BAM file and chunks the genome into smaller pieces for analysis.  Bash scripts are created in the ANALYSIS_PATH/job_scripts directory. Must be run for all diploid chromosomes (1-22 for males, 1-22 + X for females).
  * Example: 
```
		./lira split --config config-example.txt --chr 1
		
		#Run jobs in parallel for all chromosomes
		./lira split --config config-example.txt --parallel
		
		#Run jobs in parallel for all chromosomes, overwriting existing results
		./lira split --config config-example.txt --parallel --overwrite
```

3. plink (single cell and bulk)
  * Arguments: Config file
  * Runs scripts generated in (2) in parallel.
  * Example: 
```
		./lira plink --config config-example.txt
		
    #Run jobs in parallel for all chromosomes, overwriting existing results
		./lira plink --config config-example.txt --overwrite
		
		#Run jobs in parallel for all chromosomes, overwriting existing results and using a smaller batch size for jobs.
		./lira plink --config config-example.txt --batch_size 5 --overwrite
```

4. compare
  * Arguments: Single cell config file, bulk config file, chromosome (optional)
  * Compares linked variants in single cell and bulk data over a single diploid chromosome. Must be run for all diploid chromosomes.  Creates power analysis scripts run in (5).
  * Example: 
```
		./lira compare --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --chr 1
		
		#Run jobs in parallel for all chromosomes, overwriting existing results 
		./lira compare --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --parallel --overwrite
		
		#Run jobs in parallel for all chromosomes.  Wait for the completion of bulk jobs before executing.
		./lira compare --single_cell_config config-example.txt  --bulk_config config-example-bulk.txt --parallel --wait
```

5. ppower 
  * Each script analyzes the power LiRA has to detect somatic variation over different regions of the genome.
  * Example:
```
		./lira ppower --single_cell_config config-example.txt --bulk_config config-example-bulk.txt 
		
    #Run jobs in parallel for all chromosomes, overwriting existing results
		./lira ppower --single_cell_config config-example.txt --bulk_config config-example-bulk.txt --overwrite
		
		#Run jobs in parallel for all chromosomes, overwriting existing results and using a smaller batch size for jobs.
		./lira ppower --single_cell_config config-example.txt --bulk_config config-example-bulk.txt --batch_size 5 --overwrite
```

6. varcall
  * Arguments: Single cell config file, bulk config file
  * Call sSNVs and calculate genome-wide mutation rate.  Produces a VCF file `ANALYSIS_PATH/out.vcf.gz` describing the observed and passing sSNVs, as well as the germline polymorphisms used for phasing.  Also produces a rate plot `rate-plot.<bulk name>.rda` showing the sSNV rate vs. composite coverage and the threshold decision.
  * Example: 
```
    ./lira varcall --single_cell_config config-example.txt --bulk_config config-example-bulk.txt
    
    #Overwrite existing results
    ./lira varcall --single_cell_config config-example.txt --bulk_config config-example-bulk.txt --overwrite
```

7. joint (planned, not yet implemented)
  * Arguments: Single cell config file list (text file, 1 per line), bulk config file
  * Compare evidence for sSNVs across single cells analyzed against a common bulk sample and called in the same original VCF file. Produces a multi-sample VCF file `DIRECTORY/out.vcf.gz`.
  * Example: 
```
    ./lira joint --config_list config-list.txt --bulk_config config-example-bulk.txt -d joint-call
```

## Running an analysis
Suppose you have data for single cells (cells 1 through N, named sc1, sc2, ..., scN) and a corresponding bulk sample 'bulk'.  After creating config files for each single cell (sc1-N.txt) and the bulk (bulk.txt), running LiRA would amount to executing the following commands for each single cell (exemplified below for sc1):

```
    ./lira setup -c sc1.txt
    ./lira split -c sc1.txt -p
    ./lira plink -c sc1.txt
    ./lira compare -s sc1.txt -b bulk.txt -p -w
    ./lira ppower -s sc1.txt -b bulk.txt
    ./lira varcall -s sc1.txt -b bulk.txt
```

And the following commands for bulk:

```
    ./lira setup -c bulk.txt
    ./lira split -c bulk.txt -p
    ./lira plink -c bulk.txt
```

These scripts can be run simultaneously.  The single-cell `compare` step utilizing the bulk data will wait for the bulk data to complete rather than throwing an error (given the option `-w` is provided). While LiRA attempts to partition the genome in a way that balances the workload for parallel jobs, sometimes the parallel jobs for `plink` and `ppower` will not finish in the time allotted, and if this the case try rerunning with smaller `plink` and/or `ppower` batch sizes.  As LiRA will not rerun steps that have already completed, running `plink` and `ppower` with smaller batch sizes can be included in the above workflow scripts and will not run unless they need to:

Single cell:
```
    ./lira setup -c sc1.txt
    ./lira split -c sc1.txt -p
    ./lira plink -c sc1.txt
    ./lira plink -c sc1.txt -a 5 #won't run unless jobs from previous step don't finish in time
    ./lira compare -s sc1.txt -b bulk.txt -p -w
    ./lira ppower -s sc1.txt -b bulk.txt
    ./lira ppower -s sc1.txt -b bulk.txt -a 5 #won't run unless jobs from previous step don't finish in time
    ./lira varcall -s sc1.txt -b bulk.txt
```

Bulk: 
```
    ./lira setup -c bulk.txt
    ./lira split -c bulk.txt -p
    ./lira plink -c bulk.txt
    ./lira plink -c bulk.txt -a 5 #won't run unless jobs from previous step don't finish in time
```

