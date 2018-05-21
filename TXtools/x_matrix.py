import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool


def unwrap_collect_barcode_for_regions(arg, **kwarg):
    return XMatrix.collect_barcode_for_block(*arg, **kwarg)


def unwrap_cnt_barcode_for_chrm(arg, **kwarg):
    return XMatrix.cnt_barcode_for_chrm(*arg, **kwarg)


class XMatrix():
    ##comparison is between bins, by default each bin is ~10K
    def __init__(self, sf_bam, n_bins_block, bin_size, i_barcode_cutoff, sf_working_folder, sf_bed):
        self.sf_bam = sf_bam
        self.bin_size = bin_size
        self.block_size = n_bins_block * bin_size
        self.n_bins_block = n_bins_block  # # of bins per block
        self.sf_bed = sf_bed
        self.sf_working_folder = sf_working_folder
        if self.sf_working_folder[-1] != "/":
            self.sf_working_folder += "/"
        self.barcode_cutoff = i_barcode_cutoff


    def break_chrm_to_blocks(self, chrm, chrm_length):
        # break a whole chrom to blocks
        # then, break block to bins
        n_blocks = chrm_length / self.block_size
        if chrm_length % self.block_size != 0:
            n_blocks += 1
        l_blocks = []
        i_start = 0
        i_end = 0
        while i_end < chrm_length:
            i_start = i_end + 1
            i_end += self.block_size
            if i_end > chrm_length:
                i_end = chrm_length
            l_blocks.append((chrm, i_start, i_end))

        return l_blocks

    def cnt_shared_barcodes(self, n_jobs, sf_out):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        references = samfile.references
        ref_lenths = samfile.lengths
        for chrm, chrm_length in zip(references, ref_lenths):
            l_chrm_blocks = self.break_chrm_to_blocks(chrm, chrm_length)
            # for block in l_chrm_blocks:# for each block generate a file that contain the barcords of each bin (per line)
            pool = Pool(n_jobs)
            pool.map(unwrap_collect_barcode_for_regions, zip([self] * len(l_chrm_blocks), l_chrm_blocks), 1)
            pool.close()
            pool.join()
        samfile.close()

        for chrm1, chrm_length1 in zip(references, ref_lenths):
            l_chrm1_blocks = self.break_chrm_to_blocks(chrm1, chrm_length1)
            l_chrm1_with_all = []
            for chrm2, chrm_length2 in zip(references, ref_lenths):
                if chrm2 < chrm1:
                    continue
                l_chrm2_blocks = self.break_chrm_to_blocks(chrm2, chrm_length2)
                l_chrm1_with_all.append((l_chrm1_blocks, l_chrm2_blocks, chrm1, chrm2))

            pool = Pool(n_jobs)
            pool.map(unwrap_cnt_barcode_for_chrm, zip([self] * len(l_chrm1_with_all), l_chrm1_with_all), 1)
            pool.close()
            pool.join()

        #merge all the chr1_vs_chr2 files to one single file
        with open(sf_out, "w") as fout:
            for chrm1 in references:
                for chrm2 in references:
                    if chrm2 < chrm1:
                        continue
                    sf_chrm_vs_chrm_out = self.sf_working_folder + "{0}_{1}_shared_barcode_count.txt".format(chrm1,
                                                                                                             chrm2)
                    if os.path.isfile(sf_chrm_vs_chrm_out)==False:
                        print "File {0} doesn't exist!".format(sf_chrm_vs_chrm_out)
                    with open(sf_chrm_vs_chrm_out) as fin_info:
                        for line in fin_info:
                            fout.write(line)

    def cnt_barcode_for_chrm(self, record):
        l_chrm1_blocks = record[0]
        l_chrm2_blocks = record[1]
        gchrm1 = record[2]
        gchrm2 = record[3]
        sf_chrm_vs_chrm_out = self.sf_working_folder + "{0}_{1}_shared_barcode_count.txt".format(gchrm1, gchrm2)
        for chr1_block in l_chrm1_blocks:
            chr1 = chr1_block[0]
            chr1_start = chr1_block[1]
            chr1_end = chr1_block[2]
            sf_region1 = self.sf_working_folder + "{0}_{1}_{2}_block_barcodes.txt".format(chr1, chr1_start, chr1_end)
            if os.path.isfile(sf_region1) == False:
                print "File {0} doesn't exist!!!!".format(sf_region1)
                continue
            for chr2_block in l_chrm2_blocks:
                chr2 = chr2_block[0]
                chr2_start = chr2_block[1]
                chr2_end = chr2_block[2]
                sf_region2 = self.sf_working_folder + "{0}_{1}_{2}_block_barcodes.txt".format(chr2, chr2_start,
                                                                                              chr2_end)
                if os.path.isfile(sf_region2) == False:
                    print "File {0} doesn't exist!!!!".format(sf_region2)
                    continue
                self.cnt_shared_barcode_two_blocks(sf_region1, sf_region2, sf_chrm_vs_chrm_out)

    def reload_bins_from_file(self, sf_block):
        l_block_barcode = []
        with open(sf_block) as fin_block:
            for bin1 in fin_block:
                fields = bin1.split()
                chrm1 = fields[0]
                istart1 = fields[1]
                iend1 = fields[2]
                barcode_set1 = set()
                for barcode in fields[3:]:
                    barcode_set1.add(barcode)
                l_block_barcode.append((chrm1, istart1, iend1, barcode_set1))
        return l_block_barcode

    def cnt_shared_barcode_two_blocks(self, sf_block1, sf_block2, sf_out):
        with open(sf_out, "a") as fout_cnt:
            l_block1_barcode = self.reload_bins_from_file(sf_block1)
            l_block2_barcode = self.reload_bins_from_file(sf_block2)
            for record1 in l_block1_barcode:
                for record2 in l_block2_barcode:
                    bin1 = record1[3]
                    bin2 = record2[3]
                    cnt_shared = len(bin1 & bin2)
                    chrm1 = record1[0]
                    istart1 = record1[1]
                    iend1 = record1[2]
                    chrm2 = record2[0]
                    istart2 = record2[1]
                    iend2 = record2[2]
                    s_info = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(chrm1, istart1, iend1, chrm2, istart2, iend2,
                                                                          cnt_shared)
                    fout_cnt.write(s_info)
                    # each line represents a bin

    ##collect all the barcodes of a region
    def collect_barcode_for_block(self, record):
        chrm = record[0]
        istart = record[1]  # block-start
        iend = record[2]  # block-end

        l_barcode_sets = []
        for i in range(self.n_bins_block):
            set_barcode = set()
            l_barcode_sets.append(set_barcode)

        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        for alignmt in samfile.fetch(chrm, istart, iend):
            if alignmt.has_tag("BX") == False:
                continue
            pos = alignmt.pos
            offset = pos - istart
            i_bin = offset / self.bin_size
            s_barcode = alignmt.get_tag("BX")
            l_barcode_sets[i_bin].add(s_barcode)
        samfile.close()

        sf_region_out = self.sf_working_folder + "{0}_{1}_{2}_block_barcodes.txt".format(chrm, istart, iend)
        with open(sf_region_out, "w") as fout_block_barcodes:
            itemp_start = istart
            itemp_end = istart
            for barcode_set in l_barcode_sets:
                itemp_start = itemp_end + 1
                itemp_end += self.bin_size
                if len(barcode_set) > self.barcode_cutoff:
                    continue
                s_region = "{0}\t{1}\t{2}\t".format(chrm, itemp_start, itemp_end)
                fout_block_barcodes.write(s_region)
                for barcode in barcode_set:
                    fout_block_barcodes.write(barcode + "\t")
                fout_block_barcodes.write("\n")


####
if __name__ == '__main__':
    sf_bam = sys.argv[1]
    n_bins_block = 250
    bin_size = 10000
    i_barcode_cutoff = 5000
    sf_working_folder = "./tmp/"
    sf_bed = ""
    n_cores = 10
    sf_out="matrix.bed"
    xmatrix = XMatrix(sf_bam, n_bins_block, bin_size, i_barcode_cutoff, sf_working_folder, sf_bed)
    xmatrix.cnt_shared_barcodes(n_cores, sf_out)
