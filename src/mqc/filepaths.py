from mqc.utils import sel_expand
from pathlib import Path

qc_stats_dir = Path('qc_stats')

# M-bias stats processing variants
_mbias_stats_trunk = qc_stats_dir / (
    'mbias_stats/data/{samplename}_{motifs_str}'
    '_{type}_masking-{masking}')
mbias_stats = sel_expand(_mbias_stats_trunk, type='full', masking='false')
mbias_stats_masked = sel_expand(_mbias_stats_trunk, type='full', masking='true')
mbias_stats_classic = sel_expand(_mbias_stats_trunk, type='classic', masking='false')
mbias_stats_phred_threshold = sel_expand(_mbias_stats_trunk, type='phred-threshold', masking='false')
mbias_stats_phred_threshold_masked = sel_expand(_mbias_stats_trunk, type='phred-threshold', masking='true')
aggregated_mbias_stats_dir = qc_stats_dir / 'mbias_stats/aggregated_stats_for_plots/'

mbias_plots_trunk = qc_stats_dir / (
    'mbias_stats/plots/mbias-plot_'
)
qc_report_dir = qc_stats_dir / 'mbias_stats/plots'


# mbias_counts                  = "{qc_stats_dir}/{name}_mbias-counts_{motifs_str}"
# mbias_stats_p                 = "{qc_stats_dir}/{name}_mbias-stats_{motifs_str}.p"
# mbias_stats_tsv               = "{qc_stats_dir}/{name}_mbias-stats_{motifs_str}.tsv"
# mbias_stats_masked_p          = "{qc_stats_dir}/{name}_mbias-stats_masked_{motifs_str}.p"
# mbias_stats_masked_tsv        = "{qc_stats_dir}/{name}_mbias-stats_masked_{motifs_str}.tsv"
# mbias_stats_classic_trunk         = "{qc_stats_dir}/{name}_mbias-stats_{motifs_str}_classic"
# mbias_stats_classic_p         = "{qc_stats_dir}/{name}_mbias-stats_{motifs_str}_classic.p"
# mbias_stats_phred_threshold_trunk = "{qc_stats_dir}/{name}_mbias-stats_phred-threshold_{motifs_str}.p"
# mbias_stats_masked_phred_threshold_trunk = "{qc_stats_dir}/{name}_mbias-stats_masked_phred-threshold_{motifs_str}.p"
# mbias_stats_heatmap_trunk   = "{qc_stats_dir}/{name}_{motifs_str}_mbias-stats_heatmap"
# adjusted_cutting_sites_df_trunk = "{qc_stats_dir}/{name}_adjusted_cutting_sites_df_{motifs_str}"
# adjusted_cutting_sites_obj_p  = "{qc_stats_dir}/{name}_adjusted_cutting_sites_obj_{motifs_str}.p"
# adjusted_cutting_sites_df_p   = "{qc_stats_dir}/{name}_adjusted_cutting_sites_df_{motifs_str}.p"
# adjusted_cutting_sites_df_tsv = "{qc_stats_dir}/{name}_adjusted_cutting_sites_df_{motifs_str}.tsv"
# mbias_plots_trunk             = "{qc_stats_dir}/{name}_mbias-line-plot_{motifs_str}"
# cg_occurence_plot_trunk       = "{qc_stats_dir}/{name}_freq-line-plot_{motifs_str}"
# freq_with_agg_pos_plot        = "{qc_stats_dir}/{name}_freq-agg-pos-plot_{motifs_str}"
# adj_cutting_sites_plot        = "{qc_stats_dir}/{name}_adjusted-cutting-sites_barplot_{motifs_str}.png"
# cov_counts                    = "{qc_stats_dir}/{name}_coverage-counts_{motifs_str}"
# coverage_hist                 = "{qc_stats_dir}/{name}_coverage-hist_{motifs_str}.png"
# meth_calls_basepath           = "meth_calls/mcalls"
