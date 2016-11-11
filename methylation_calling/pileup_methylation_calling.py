import mqc
import itertools


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
