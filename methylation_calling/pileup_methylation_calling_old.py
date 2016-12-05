import mqc
import mqc.methylation_calling.segment_methylation
from mqc.conv_err import call_frag_conversion_error


class MotifMethylationCaller:
    def __init__(self,
                 max_number_of_unconv_cytosines,
                 trimming_site_array=None,
                 cutting_sites_dict=None,
                 ):
        self.max_number_of_unconv_cytosines = max_number_of_unconv_cytosines
        self.trimming_site_array = trimming_site_array
        self.cutting_sites_dict = cutting_sites_dict

    def get_total_methylation_stats_at_pileup(self, motif_pileups,
                                              index_position: 'mqc.IndexPosition'):
        meth_counts_dict = {'methylated': 0,
                            'unmethylated': 0}

        watson_motif_seq = index_position.watson_motif

        for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
            if motif_base in ['C', 'G']:
                self.get_meth_counts_at_pos(pileup_reads,
                                            watson_ref_base=motif_base,
                                            index_position=index_position)

        try:
            beta = meth_counts_dict['methylated'] / (
                meth_counts_dict['methylated'] + meth_counts_dict['unmethylated'])
        except ZeroDivisionError:
            beta = -1

        return beta, meth_counts_dict['methylated'], meth_counts_dict['unmethylated']

    def get_meth_counts_at_pos(self, pileup_reads, watson_ref_base,
                               index_position: 'mqc.IndexPosition',
                               meth_counts_dict):

        sorted_pileup_reads = sorted(pileup_reads, key=lambda read: read.alignment.query_name)

        read_in_catch = False
        for pileup_read in sorted_pileup_reads:
            if not pileup_read.query_position:
                if read_in_catch:
                    add_meth_status_of_non_overlapping_read(cached_read, meth)
                    read_in_catch = False
            elif call_frag_conversion_error(pileup_read, index_position,
                                            self.max_number_of_unconv_cytosines):
                if read_in_catch:
                    add_meth_status_of_non_overlapping_read(cached_read)
                    read_in_catch = False
            elif call_whether_has_to_be_trimmed(pileup_read,
                                                trimming_site_array=self.trimming_site_array):
                if read_in_catch:
                    add_meth_status_of_non_overlapping_read(cached_read, meth_counts_dict,
                                                            meth_counts_dict=meth_counts_dict)
                    read_in_catch = False
            elif read_in_catch:
                if cached_name == pileup_read.alignment.query_name:
                    call_meth_from_overlapping_reads(cached_read, pileup_read,
                                                     watson_ref_base=watson_ref_base,
                                                     meth_counts_dict=meth_counts_dict)
                    read_in_catch = False
                else:
                    add_meth_status_of_non_overlapping_read(cached_read,
                                                            watson_ref_base=watson_ref_base,
                                                            meth_counts_dict=meth_counts_dict)
                    cached_read = pileup_read
                    cached_name = pileup_read.alignment.query_name
            else:
                read_in_catch = True
                cached_name = pileup_read.alignment.query_name
                cached_read = pileup_read

        # Take care of last read in pileup
        if read_in_catch:
            add_meth_status_of_non_overlapping_read(cached_read, watson_ref_base,
                                                    meth_counts_dict=meth_counts_dict)

def add_meth_status_of_non_overlapping_read(read, meth_counts_dict):
    try:
        meth_counts_dict[read.meth_status_at_pileup_pos] += 1
    except KeyError:
        pass


def call_meth_from_overlapping_reads(cached_read, pileup_read, watson_ref_base, meth_counts_dict):
    # TODO: does this logic work for undirectional protocols?
    meth_status_str_1 = call_meth_at_base_except_when_N_was_observed(cached_read, watson_ref_base)
    meth_status_str_2 = call_meth_at_base_except_when_N_was_observed(pileup_read, watson_ref_base)
    if meth_status_str_1 == meth_status_str_2:
        try:
            meth_counts_dict[meth_status_str_1] += 1
        except KeyError:
            pass
            # else:
            #     discard fragment
