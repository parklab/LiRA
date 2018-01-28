#! /home/clb36/parkhome/installed/local/bin/python2.7
#samtools_view.py     Bernardo  3mar2016  v. 5mar  5PM
#requires python2.7  
#usage: samtools_view.py   --bam reads.bam  --fasta reference_genomic2.fasta  --chr ctg2_left_rc --start 177718 --end 177722
#output:
#ctg2_left_rc    177718  A       ctg2_100x_PB_L_3933_1_0 11659   A
#ctg2_left_rc    177718  A       ctg2_100x_PB_L_5164_1_0 None    *
#code based on: http://pysam.readthedocs.org/en/latest/
#argparse info: http://www.cyberciti.biz/faq/python-command-line-arguments-argv-example/
import pysam
import argparse
import csv
parser = argparse.ArgumentParser(description='usage: samtools_view.py --bam reads.bam --bed bed_file.bed')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bedfile name',required=True)
args = parser.parse_args()
bamfile = pysam.AlignmentFile(args.bam, "rb")
with open(args.bed) as bed:
	reader = csv.reader(bed, delimiter="\t")
	sites = list(reader)

for site in sites:
	pile={'A': 0, 'C': 0, 'G': 0, 'T': 0}
	start = int(site[1])
	end = int(site[2])
	pileup = bamfile.pileup(site[0],start,end,stepper="all",max_depth=500000)
	for pileupColumn in pileup:
		for pileupRead in pileupColumn.pileups:
			if (pileupColumn.pos >= start) and (pileupColumn.pos < end) and (not pileupRead.is_refskip) and (not pileupRead.is_del):
				base = pileupRead.alignment.query_sequence[pileupRead.query_position]
				if((base != 'N') and (base != '*')):
					pile[base]=pile[base] + 1
	print ('%s\t%s\t%s\t%s' % (pile['A'],pile['C'],pile['G'],pile['T']))
bamfile.close()
