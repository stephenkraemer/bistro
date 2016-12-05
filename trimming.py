def call_to_trimm_flen_specific_cutting_sites(read, trimming_site_array):
    try:
        tlen = read.alignment.template_length
        lowest_allowed_pos = trimming_site_array[0][tlen]
        highest_allowed_pos = trimming_site_array[1][tlen]
    except IndexError:  # fragment exceptionally long
        lowest_allowed_pos = trimming_site_array[0][-1]
        highest_allowed_pos = trimming_site_array[1][-1]
    return not (lowest_allowed_pos <= read.query_position <= highest_allowed_pos)


def call_to_trimm_minimal_cutting_sites_only(read: 'pysam.PileupRead', cutting_sites_dict):
    # TODO: Improve template length determination
    pos_in_read = read.query_position
    # TODO
    # bs_seq_strand = read.get_bsseq_strand()
    bsseq_strand = 'W-BC'
    if pos_in_read < cutting_sites_dict[bsseq_strand][0]:
        return True
    elif pos_in_read > read.alignment.template_length - cutting_sites_dict[bsseq_strand][1]:
        return True


