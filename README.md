# LiRA [![Build Status](https://travis-ci.com/parklab/LiRA.svg?token=EkzyvwdZ2jcY78ErmS88&branch=master)](https://travis-ci.com/parklab/LiRA) [![codecov](https://codecov.io/gh/parklab/LiRA/branch/master/graph/badge.svg?token=PMpspJNu0Z)](https://codecov.io/gh/parklab/LiRA)

LiRA is a tool to detect single-cell specific SNVs (scsSNVs).
Its functionalities include:


# Pre-reqs

LiRA is written entirely in python and has been developed using the following python version:
- **python**: 2.7.6 (default, Aug  3 2015, 17:43:52)  [GCC 4.4.7 20120313 (Red Hat 4.4.7-11)] 

We use pip and virtual enviornments to take care of our dependency management
- [**pip**](https://pip.pypa.io/en/stable/installing/)
- [**virtualenvwrapper**](https://virtualenvwrapper.readthedocs.io/en/latest/install.html#installation) (optional, but recommended)


# Installation

`$ git clone https://github.com/parklab/LiRA.git && cd LiRA`

`$ mkvirtualenv LiRA-env` (optional, but recommended)

`$ pip install -r requirements.txt`

# Running Tests
`$ python tests.py`

```sh
$ python scripts/get_reference_set_from_fasta.py
$ ./scripts/sort_reference_sets.sh
```

These scripts will generate one file per chromosome containing the reference sets of MS repeats (coordinate-sorted).

Make sure to add the current directory to your PYTHONPATH:
export PYTHONPATH=$PYTHONPATH:$PWD

<!-- ### Detection of microsatellites

![examples_detection](https://user-images.githubusercontent.com/5588266/28093196-42fedc6a-6664-11e7-9e8a-04d555a88e7e.png)
-->

## LiRA parameters

Once the reference sets are ready, LiRA can be used.

Information on the parameters can be accessed by typing:
```sh
$ python msi_profiler.py --help
```

Example of usage:
```sh
python path_to_msi_profiler/msi_profiler.py  --tumor_bam tumor.bam  --normal_bam normal.bam --bed hets_SNPs_chr22.bed --chromosomes 22  --fasta path_to_LiRA/chrs_fa/ --output_prefix example_chr22  --mode unphased --nprocs 8  --reference_set path_to_reference_sets_folder --min_coverage 8 --min_MS_length 6 --flank_size 10 --rus 1 2 3 4 5 6  --tolerated_mismatches 0
```


### Example of usage
```sh
python msi_profiler.py --tumor_bam test_tumor.bam  --normal_bam test_normal.bam  --bed germline_het_SNPs.bed --chromosomes 21 22 X --fasta ./chrs_fa/  --output_prefix example_unphased --mode unphased --nprocs 2  --reference_set path_to_reference_sets_folder  --min_coverage 8 --min_MS_length 6 --flank_size 5 --rus 1 2 3 4 5 6 
```

### Example of output

# Notes
Currently, LiRA considers by default the hg19 assembly of the human genome.
However, LiRA can be used with other assemblies if the bam files and the fasta sequences from which the reference sets are derived are concordant.

<!--
# How to cite
The details of LiRA have been published in:
XXX
-->

# Contact
If you have any questions or suggestions please contact us: 
- Isidro Cortes Ciriano: isidrolauscher at gmail.com  or isidro\_cortesciriano at hms.harvard.edu
- Peter J Park: peter\_park at hms.harvard.edu



