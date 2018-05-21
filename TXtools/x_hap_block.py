##02/28/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool


class Region():
    def __init__(self, chrm, start, end):  # for a whole chrm, set the "start" and "end" to be "None"
        self.chrm = chrm
        self.start = start
        self.end = end

    def getRegion(self):
        return (self.chrm, self.start, self.end)


def unwrap_self_block_statistic_by_region(arg, **kwarg):
    return XHapBlock.run_get_block_info(*arg, **kwarg)


class XHapBlock():
    def __init__(self, sf_bam, working_folder, n):
        self.sf_bam = sf_bam
        self.cores = n
        self.working_folder = working_folder
        if os.path.exists(working_folder) == False:
            print "Error: The provided temp folder doesn't exist!!!"
            self.working_folder = "./"
        elif self.working_folder[-1] != "/":
            self.working_folder += "/"

    def run_get_block_info(self, record):
        region = record[0]
        sf_bam = record[1]
        sf_out = record[2]
        chrm, istart, iend = region.getRegion()

        samfile = pysam.AlignmentFile(sf_bam, "rb")
        m_block_info = {}
        m_hap_block = {}  ###save the block start and end info
        for alignmt in samfile.fetch(chrm, istart, iend):

            pos = alignmt.pos

            i_molecule_id = -1
            if alignmt.has_tag("MI"):
                i_molecule_id = alignmt.get_tag("MI")

            s_barcode = "empty"
            if alignmt.has_tag("BX"):  # if there is a barcode
                s_barcode = alignmt.get_tag("BX")

            i_phase_block = -1
            if alignmt.has_tag("PS"):
                i_phase_block = alignmt.get_tag('PS')

            if i_phase_block != -1:
                if i_phase_block not in m_hap_block:
                    m_hap_block[i_phase_block] = []
                    m_hap_block[i_phase_block].append(pos)
                    m_hap_block[i_phase_block].append(pos)
                else:
                    tmp_start = m_hap_block[i_phase_block][0]
                    tmp_end = m_hap_block[i_phase_block][1]
                    if pos > tmp_end:
                        m_hap_block[i_phase_block][1] = pos
                    if pos < tmp_start:
                        m_hap_block[i_phase_block][0] = pos

            i_hap = -1
            if alignmt.has_tag("HP"):
                i_hap = alignmt.get_tag('HP')

            if i_phase_block == -1:
                continue

            if i_phase_block not in m_block_info:
                m_block_info[i_phase_block] = {}

            if i_molecule_id not in m_block_info[i_phase_block]:
                m_block_info[i_phase_block][i_molecule_id] = {}

            if s_barcode not in m_block_info[i_phase_block][i_molecule_id]:
                m_block_info[i_phase_block][i_molecule_id][s_barcode] = 1
            else:
                m_block_info[i_phase_block][i_molecule_id][s_barcode] += 1
        samfile.close()

        with open(sf_out, "w") as fout_block:
            for i_phase_block in m_block_info:
                for i_molecule_id in m_block_info[i_phase_block]:
                    for s_barcode in m_block_info[i_phase_block][i_molecule_id]:
                        n_cnt_read = m_block_info[i_phase_block][i_molecule_id][s_barcode]
                        istart_tmp = m_hap_block[i_phase_block][0]
                        iend_tmp = m_hap_block[i_phase_block][1]
                        s_block = "{0}\t{1}\t{2}\t".format(chrm, istart_tmp, iend_tmp)
                        s_info = "{0}\t{1}\t{2}\t{3}\n".format(i_phase_block, i_molecule_id, s_barcode, n_cnt_read)
                        fout_block.write(s_block + s_info)


    def get_block_info(self, sf_out):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        references = samfile.references
        l_chrm_records = []
        for chrm in references:
            region = Region(chrm, None, None)
            sf_region_out = self.working_folder + "{0}.block_stat.txt".format(chrm)
            l_chrm_records.append((region, self.sf_bam, sf_region_out))
        samfile.close()

        pool = Pool(self.cores)
        pool.map(unwrap_self_block_statistic_by_region, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        #merge the files:
        with open(sf_out, "w") as fout_rslt:
            fout_rslt.write("Chrom\tBlock_start\tBlock_end\tPhased_Block_Id\tGlobal_Molecule_Id\tBarcode\t#_of_reads\n")
            for chrm in references:
                sf_region_out = self.working_folder + "{0}.block_stat.txt".format(chrm)
                print sf_region_out
                if os.path.isfile(sf_region_out)==False:
                    continue
                with open(sf_region_out) as fin_region:
                    for line in fin_region:
                        fout_rslt.write(line)

