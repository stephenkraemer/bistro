import mqc
import mqc.methylation_calling.segment_methylation
from mqc.conv_err import call_frag_conversion_error

bs_strand_dict_directional_protocol = {1: {'forward': 'W-BC',
                                           'reverse': 'C-BC'},
                                       2: {'forward': 'C-BC-Rv',
                                           'reverse': 'W-BC-Rv'}}


def get_total_methylation_stats_at_pileup(motif_pileups,
                                          index_position: 'mqc.IndexPosition',
                                          max_number_of_unconv_cytosines,
                                          trimming_site_array):
    meth_counts_dict = {'methylated': 0,
                        'unmethylated': 0}

    watson_motif_seq = index_position.watson_motif

    for motif_base, pileup_reads in zip(watson_motif_seq, motif_pileups):
        if motif_base in ['C', 'G']:
            get_meth_counts_at_pos(
                pileup_reads,
                watson_ref_base=motif_base,
                index_position=index_position,
                max_number_of_unconv_cytosines=max_number_of_unconv_cytosines,
                trimming_site_array=trimming_site_array, meth_counts_dict=meth_counts_dict)

    try:
        beta = meth_counts_dict['methylated'] / (
            meth_counts_dict['methylated'] + meth_counts_dict['unmethylated'])
    except ZeroDivisionError:
        beta = -1

    return beta, meth_counts_dict['methylated'], meth_counts_dict['unmethylated']


def get_meth_counts_at_pos(pileup_reads, watson_ref_base,
                           index_position: 'mqc.IndexPosition',
                           max_number_of_unconv_cytosines,
                           trimming_site_array,
                           meth_counts_dict):

    sorted_pileup_reads = sorted(pileup_reads, key=lambda read: read.alignment.query_name)

    read_in_catch = False
    for pileup_read in sorted_pileup_reads:
        if not pileup_read.query_position:
            if read_in_catch:
                call_meth_on_non_overlapping_read(cached_read,
                                                  watson_ref_base=watson_ref_base,
                                                  meth_counts_dict=meth_counts_dict)
                read_in_catch = False
        elif call_frag_conversion_error(pileup_read, index_position, max_number_of_unconv_cytosines):
            if read_in_catch:
                call_meth_on_non_overlapping_read(cached_read,
                                                  watson_ref_base=watson_ref_base,
                                                  meth_counts_dict=meth_counts_dict)
                read_in_catch = False
        elif call_whether_has_to_be_trimmed(pileup_read,
                                            trimming_site_array=trimming_site_array):
            if read_in_catch:
                call_meth_on_non_overlapping_read(cached_read, meth_counts_dict,
                                                  meth_counts_dict=meth_counts_dict)
                read_in_catch = False
        elif read_in_catch:
            if cached_name == pileup_read.alignment.query_name:
                call_meth_from_overlapping_reads(cached_read, pileup_read,
                                                 watson_ref_base=watson_ref_base,
                                                 meth_counts_dict=meth_counts_dict)
                read_in_catch = False
            else:
                call_meth_on_non_overlapping_read(cached_read,
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
        call_meth_on_non_overlapping_read(cached_read, watson_ref_base,
                                          meth_counts_dict=meth_counts_dict)


def call_meth_on_non_overlapping_read(read, watson_ref_base, meth_counts_dict):
    meth_status_str = call_meth_at_base_except_when_N_was_observed(
        read, watson_ref_base)
    try:
        meth_counts_dict[meth_status_str] += 1
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


def call_meth_at_base_except_when_N_was_observed(pileup_read, watson_ref_base):
    direction_str = 'reverse' if pileup_read.alignment.is_reverse else 'forward'
    mate = 1 if pileup_read.alignment.is_read1 else 2
    bs_seq_strand = bs_strand_dict_directional_protocol[mate][direction_str]
    pos_in_read = pileup_read.query_position
    if not pos_in_read:
        meth_status_str = 'NA'
    else:
        observed_watson_base = pileup_read.alignment.query_sequence[pos_in_read]
        if observed_watson_base == 'N':
            meth_status_str = 'NA'
        else:
            meth_status_str = mqc.methylation_calling.segment_methylation.meth_dict[
                bs_seq_strand][watson_ref_base][observed_watson_base]
    return meth_status_str


def call_whether_has_to_be_trimmed(read, trimming_site_array):
    try:
        tlen = read.alignment.template_length
        lowest_allowed_pos = trimming_site_array[0][tlen]
        highest_allowed_pos = trimming_site_array[1][tlen]
    except IndexError:  # fragment exceptionally long
        lowest_allowed_pos = trimming_site_array[0][-1]
        highest_allowed_pos = trimming_site_array[1][-1]
    return not (lowest_allowed_pos <= read.query_position <= highest_allowed_pos)


"""
def find_next_good_read(sorted_pileup_reads, index_position,
                        max_number_of_unconv_control_cyts, trimming_site_array):
    while 1:
        cached_read = next(sorted_pileup_reads)
        has_conversion_error = call_frag_conversion_error(
            cached_read, index_position, max_number_of_unconv_control_cyts)
        has_to_be_trimmed = call_whether_has_to_be_trimmed(cached_read, trimming_site_array)
        if not (has_conversion_error or has_to_be_trimmed):
            return cached_read, cached_read.alignment.query_name



# ?
            for pileup_read in pileup_reads:
                # TODO: more robust skipping of bad reads
                if meth_status_str == 'methylated':
                    n_meth += 1
                elif meth_status_str == 'unmethylated':
                    n_unmeth += 1


def get_total_meth_stats_at_pileup(meth_pileup: 'mqc.MethPileup', trimming_mode):
    n_meth, n_unmeth = 0, 0
    for pileup_segment in itertools.chain.from_iterable(meth_pileup.pileup_segments_per_pos):
        if not pileup_segment.is_ok:
            continue

        if pileup_segment.omit_due_to_overlap:
            continue

        if pileup_segment.omit_due_to_trimming_dict[trimming_mode]:
            continue

        if pileup_segment.has_frag_conversion_error:
            continue

        if pileup_segment.meth_status_str == 'methylated':
            n_meth += 1
        elif pileup_segment.meth_status_str == 'unmethylated':
            n_unmeth += 1

    n_total = n_meth + n_unmeth
    try:
        beta = n_meth / n_total
    except ZeroDivisionError:
        beta = None

    return beta, n_meth, n_total
"""
