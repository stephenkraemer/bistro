import mqc


def call():
    index_file = mqc.IndexFile('/home/stephen/projects/mqc/mqc/test/data/'
                               'chr11_16815793-16824254.cg.bed.edited.gz')
    pileups = mqc.PysamPileupEngine(index_file,
                                    bam_abspath='/home/stephen/projects/mqc/mqc/test/data/'
                                                'b_cells_rep1_chr11_16815793-16824254.bam',
                                    mapq_threshold=1,
                                    phred_score_threshold=1,
                                    run_mode='call',
                                    overlap_handling_mode='discard_disagreeing_overlaps')
    for pileup in pileups:
        print(pileup.index_position.chrom,
              pileup.index_position.start)
        # print(pileup)
        beta, n_meth, n_unmeth = pileup.get_total_meth_stats()
        print(beta, n_meth, n_unmeth)

        try:
            c += 1
        except UnboundLocalError:
            c = 1
        if c > 2:
            break

call()
