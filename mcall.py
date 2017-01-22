def tabulate_meth_calls(bam_path, index_file_path, output_file, cutting_site_array,
                        max_flen_considered_for_trimming):
    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)

    with gzip.open(output_file, 'wt') as fout:
        for motif_pileups, curr_idx_pos in motif_pileup_iter:
            beta, n_meth, n_unmeth = mc.call_meth_at_pileup(
                    motif_pileups,
                    index_position=curr_idx_pos,
                    cutting_site_array=cutting_site_array,
                    max_flen_considered_for_trimming=max_flen_considered_for_trimming)

            print(curr_idx_pos.chrom,
                  curr_idx_pos.start,
                  curr_idx_pos.end,
                  beta,
                  n_meth,
                  n_unmeth,
                  sep='\t', end='\n', file=fout)
