import mqc.trimming

import mqc
import mqc.qc_taggers.overlap


def call_meth_at_pileup(motif_pileups, index_position: 'mqc.IndexPosition', cutting_site_array,
                        max_flen_considered_for_trimming):
    n_meth = 0
    n_unmeth = 0

    watson_motif_seq = index_position.watson_motif
    for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
        if motif_base in ['C', 'G']:
            mqc.qc_taggers.trimming.set_trimming_flag(pileup_reads, cutting_site_array,
                                                      max_flen_considered_for_trimming=max_flen_considered_for_trimming)
            mqc.qc_taggers.overlap.tag_overlaps(pileup_reads)
            for read in pileup_reads:
                # read: mqc.BSSeqPileupRead
                if read.qc_fail_flag:
                    continue
                if read.overlap_flag:
                    continue
                if read.trimm_flag:
                    continue
                meth_status_flag = read.get_meth_status_at_pileup_pos(motif_base)
                if meth_status_flag == 8:
                    n_meth += 1
                elif meth_status_flag == 4:
                    n_unmeth += 1

    try:
        beta = n_meth / (n_meth + n_unmeth)
    except ZeroDivisionError:
        beta = 'NA'

    return beta, n_meth, n_unmeth
