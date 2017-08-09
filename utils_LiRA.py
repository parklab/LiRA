import bisect
import csv

import pysam
import numpy as np

def phased(msi_obj, sites, bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        dict_out = {}
        for site in sites:
            start = int(site[1])
            end = int(site[2])
            chr = str(site[0])
            bases = [site[3], site[4]]
            reads = [
                read for read in bamfile.fetch(
                    chr,
                    start,
                    end,
                    multiple_iterators=True
                )
                ]  ## keep this as is, do not put conditions inside here

            reads = [
                read for read in reads if read.is_proper_pair and
                read.is_duplicate == False and
                read.mapping_quality >= msi_obj.mapping_quality
                ]
            mates = [
                bamfile.mate(read) for read in reads 
                ]
			reads = reads + mates
            if len(reads) > msi_obj.min_coverage:
                for read in reads:
                    read_sequence = read.seq
                    reps = find_repeats(
                        read_sequence,
                        msi_obj.flank_size,
                        msi_obj.repeat_units
                    )
                    if len(reps) > 0:
                        # get the SNP allele in this read
                        start_read = read.reference_start
                        end_read = read.reference_end
                        aligned_pos = read.get_reference_positions(
                            full_length=True
                        )  # True) reports none for  soft-clipped positions
                        try:
                            idx = aligned_pos.index(start)
                        except:
                            continue
                        snp_read = read_sequence[idx]
                        if snp_read not in bases:
                            continue
                        for microsatellite in reps:
                            rs = microsatellite[1]
                            re = microsatellite[2]
                            difference = re - rs + 1
                            # use the reference set here to get the
                            # position on the right
                            ini = start_read + rs  #
                            idx2 = binary_search(
                                msi_obj.reference_set_ini_end_dict[chr],
                                str(ini + 1)
                            )

                            if idx2 == -1:
                                continue
                            refset_now = msi_obj.reference_set_dict[chr][idx2]
                            diff_ref = int(refset_now[2]) - int(
                                refset_now[1]) + 1
                            with pysam.FastaFile(
                                    filename=msi_obj.fasta_dict[chr]
                            ) as fasta_file:
                                flank_right_ref = fasta_file.fetch(
                                    "chr" + str(site[0]), ini + diff_ref,
                                    ini + diff_ref + msi_obj.flank_size).upper()
                                flank_left_ref = fasta_file.fetch(
                                    "chr" + str(site[0]), ini -
                                    msi_obj.flank_size,
                                    ini).upper()
                            posfl = (start_read + rs - msi_obj.flank_size)
                            if posfl >= start_read:
                                flank_left = read_sequence[
                                             rs - msi_obj.flank_size:rs]
                                mismatches_left = sum(a != b for a, b in
                                                      zip(flank_left,
                                                          flank_left_ref))
                            else:
                                mismatches_left = 10000
                            posflr = start_read + re + msi_obj.flank_size
                            if posflr <= end_read:
                                flank_right = read_sequence[
                                              re + 1:re + 1 +
                                                     msi_obj.flank_size]
                                mismatches_right = sum(a != b for a, b in
                                                       zip(flank_right,
                                                           flank_right_ref))
                            else:
                                mismatches_right = 10000
                            mismatches = mismatches_left + mismatches_right
                            if mismatches <= msi_obj.tolerated_mismatches:
                                key_now = site[0] + "\t" + str(ini) + "\t" + \
                                          refset_now[3] + "\t" + refset_now[
                                              4] + "\t" + refset_now[
                                              5] + "\t" + refset_now[
                                              6] + "\t" + snp_read + "\t" + \
                                          str(site[1])
                                if dict_out.has_key(key_now):
                                    dict_out[key_now] = np.append(
                                        dict_out[key_now], difference)
                                else:
                                    dict_out[key_now] = difference

    return dict_out


def multiprocessing_lock_init(l):
    global lock
    lock = l
