""" Default python config file

Why a python config file?
-------------------------
I think JSON and TOML are the main alternatives
which come to mind. I actually tried TOML and aborted due to parser issues
with both toml and pytoml packages.
  - JSON: no comments, no trailing commas. Comments are crucial for
  development and documentation of the plotting spec. We could use
  minify, but that seems convoluted
  - TOML: no stable parser for nested inline dicts at the time of writing,
          the standard way of writing dicts requires abstraction from
          the python dict we are actually trying to construct here
  - YAML: too complicated
  - XML: too complicated
y_breaks: Union[List[int], int]
x_breaks: Union[List[int], int]
y_lim: 'auto', or tuple of ints

add option to aggregate levels? eg. aggregate seq_contexts?

groupby several levels


allow window size for ylimits, e.g. beta_values= 0.3 to place the 0.3 window
which contains most points? or to place window of size max(0.3, size needed to accomodate all points)

never share y limits across different variables

Implementation checklist
- breaks can be integer or list of int
- axis limits can be tuple or string (and in future perhaps window size as float)
Format
======

Caveats
-------
- don't use tuples to group datasets (use lists)


      
      
"""

default_config = {

    # Defaults
    # ==================================================================
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
                       'share': True,
                       'rotate_labels': False,
                       },
            'x_axis': {'breaks': [0, 30, 60, 90, 120, 150],
                       'limits': {
                           'phred': 45,
                       },
                       'share': True,
                       'rotate_labels': True},
            'plot': ['line'],
        },

        "post_agg_filters": {
            'flen': [60, 75, 90, 100, 120, 150, 190, 290, 490],
            'phred': [9, 14, 19, 24, 29, 34, 39, 40],
        },

    },

    # Plot groups
    # ==================================================================


    # Positional bias plots
    # ------------------------------------------------------------------

    "Positional_bias_plots": {
        "datasets": [["full", "trimmed"]],
        "pre_agg_filters": {
            'motif': ['CG']
        },
        "aes_mappings": [

            # Global M-bias plot before and after trimming
            {
                "x": "pos",
                "y": "beta_value",
                "row": "dataset",
                "column": "bs_strand",
                "color": None
            },

            # Fragment length-stratified M-bias plot
            {
                "x": "pos",
                "y": "beta_value",
                "row": "dataset",
                "column": "bs_strand",
                "color": "flen"
            },

            {
                "x": "pos",
                "y": "beta_value",
                "row": "dataset",
                "column": "bs_strand",
                "color": "phred"
            },

            {
                "x": "pos",
                "y": "beta_value",
                "row": "dataset",
                "column": "bs_strand",
                "color": "seq_context"
            },
        ],
    },

    # Mate level positional bias plots
    # ------------------------------------------------------------------
    # "mate-level plots": {
    #     "datasets": [["full", "trimmed"]],
    #     "pre_agg_filters": {
    #         'motif': ['CG']
    #     },
    #     "aes_mappings": [],
    # },

    # # Coverage plots
    # # ------------------------------------------------------------------
    # 'Coverage plots': {
    #     "datasets": [["full", "trimmed"]],
    #     "pre_agg_filters": {
    #         'motif': ['CG']
    #     },
    #     "plot_params": {
    #         'sharey': False,
    #     },
    #     "aes_mappings": [
    #
    #         # M-bias plot
    #         {
    #             'x': 'pos',
    #             'y': 'value',
    #             'color': None,
    #             "row": 'dataset',
    #             "column": "bs_strand"
    #         },
    #
    #         # Fragment-length stratified M-bias plot
    #         {
    #             'x': 'pos',
    #             'y': 'value',
    #             'color': 'flen',
    #             "row": 'dataset',
    #             "column": "bs_strand"
    #         },
    #
    #     ],
    # },

    # # BS-strand or mate methylation bias
    # # ------------------------------------------------------------------
    # # See paper about mate bias due to software version
    # "strand_mate_meth_bias": {
    #     "datasets": ["full_phred-threshold",
    #                  "trimmed_phred-threshold"],
    #     "aes_mappings": [
    #
    #         {
    #             'x': 'bs_strand',
    #             'y': 'beta_value',
    #             # 'color': [None, 'phred'],
    #             'color': 'phred',
    #             'row': 'dataset',
    #         },
    #
    #         {
    #             'x': 'bs_strand',
    #             'y': 'beta_value',
    #             'color': 'phred',
    #             'row': 'dataset',
    #             'column': 'seq_context',
    #         },
    #
    #         # {
    #         #     'x': 'bs_strand',
    #         #     'y': 'value',
    #         #     'color': 'phred',
    #         #     'row': 'dataset',
    #         #     'column': 'statistic',
    #         # },
    #         #
    #         #
    #         #
    #
    #     ],
    # },

    # # Sequence context bias due to phred filtering
    # # ------------------------------------------------------------------
    # # - for comparison of general sequence context bias, one could also
    # #   compare multiple samples
    #
    # # TODO: add phred-trheshold dfs in demo snakefile
    # "seq_context_plots": {
    #     'plot_params': {
    #         'plot': ['bar', ['line', 'point']],
    #         'x_axis': {'breaks': 'auto'}
    #     },
    #     "datasets": ["full_phred-threshold",
    #                  "trimmed_phred-threshold"],
    #     "aes_mappings": [
    #
    #         {
    #             "x": "seq_context",
    #             "y": "beta_value",
    #             "row": "dataset",
    #             "column": "bs_strand",
    #             "color": "phred"
    #         },
    #
    #         # {
    #         #     "x": "seq_context",
    #         #     "y": "values",
    #         #     "row": 'bs_strand',
    #         #     "column": 'statistic',
    #         #     "color": "phred"
    #         # },
    #         #
    #     ],
    # }

}
