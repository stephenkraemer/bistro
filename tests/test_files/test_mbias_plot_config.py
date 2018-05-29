default_config = {
    "defaults": {
        'plot_params': {
            "panel_height_cm": 6,
            "panel_width_cm": 6,
            "theme": "paper",
            'y_axis': {'breaks': 5,
                       'limits': {
                           'beta_value': [(0, 1), (0.5, 1), 'auto'],
                           'n_meth': 'auto',
                           'n_unmeth': 'auto', },
                       'share': True
                       },
            'x_axis': {'breaks': [0, 30, 60, 90, 120, 150],
                       'share': True,
                       'rotate_labels': True},
            'plot': ['line'],
        },
        "pre_plot_filters": {
            'flen': [60, 75, 90, 100, 120, 150, 190, 290, 490],
            'phred': [0, 5, 10, 15, 20, 25, 30, 35, 40],
        },
    },
    "group1": {
        "datasets": [["full", "trimmed"]],
        # necessary because currently only trimming for CG
        "pre_agg_filters": {
            'motif': ['CG']
        },
        "aes_mappings": [
            {'x': 'pos', 'y': 'beta_value',
             "row": "dataset", "column": "bs_strand", "color": None},
            {'x': 'pos', 'y': 'beta_value',
             "row": "dataset", "column": "bs_strand", "color": "flen"},
            {'x': 'pos', 'y': 'beta_value',
             "row": "dataset", "column": "bs_strand", "color": None},
        ],
    },
    "group2": {
        "datasets": ["full", "trimmed"],
        # necessary because currently only trimming for CG
        "aes_mappings": [
            {'x': 'pos', 'y': 'beta_value',
             "row": "motif", "column": "bs_strand", "color": None},
        ],
    },
}
