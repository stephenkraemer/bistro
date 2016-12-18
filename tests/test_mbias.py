import pytoml
import pickle
import mqc.mbias

load_old_mbias_data = True
load_cutting_sites = True

with open('./config.default.toml') as f_toml:
    config_dict = pytoml.load(f_toml)

mbias_stats_array_p = '/home/kraemers/projects/mqc/mqc/test/results/mbias_stats_array.p'
with open(mbias_stats_array_p, 'rb') as f:
    mbias_stats_arr = pickle.load(f)

minimal_cutting_sites = mqc.mbias.MinimalCuttingSites(
    relative_cutting_site_dict=config_dict['trimming']['relative_to_fragment_ends'],
    max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming'],
    max_read_length_bp=config_dict['data_properties']['max_read_length_bp'])

if load_old_mbias_data:
    with open('/home/kraemers/projects/mqc/mqc/test/results/mbias_data.p', 'rb') as f:
        mbias_data = pickle.load(f)
else:
    mbias_data = mqc.mbias.MbiasData(
        max_flen_considered_for_trimming=config_dict['trimming']['max_flen_considered_for_trimming'],
        max_read_length_bp=config_dict['data_properties']['max_read_length_bp'],
        required_n_events_for_cutting_site_determination=(
            config_dict['trimming']['required_n_events_for_cutting_site_determination']),
        max_window_size_for_smoothing=config_dict['trimming']['max_window_size_for_smoothing'],
        mbias_stats_array=mbias_stats_arr)
    # Not Implemented:
    # minimal_cutting_sites_df=minimal_cutting_sites.get_df())
    with open('/home/kraemers/projects/mqc/mqc/test/results/mbias_data.p', 'wb') as f:
        pickle.dump(mbias_data, f)

if load_cutting_sites:
    with open('/home/kraemers/projects/mqc/mqc/test/results/cutting_sites.p', 'rb') as f:
        cutting_sites = pickle.load(f)
        cutting_sites_df = cutting_sites.get_df()
else:
    cutting_sites = mqc.mbias.AdjustedMbiasCuttingSites(mbias_data, calling_mode='standard',
                                                        min_len=50, max_std=0.1)
    cutting_sites_df = cutting_sites.get_df()
    with open('/home/kraemers/projects/mqc/mqc/test/results/cutting_sites.p', 'wb') as f:
        pickle.dump(cutting_sites, f)

mqc.mbias.cutting_sites_area_plot(mbias_cutting_sites_df=cutting_sites_df,
                                  output_path='/home/kraemers/temp/cutting_sites.png')

mbias_data._adjusted_cutting_sites_df = cutting_sites.get_df()
plotter = mqc.mbias.MbiasDataPlotter(mbias_data=mbias_data)

plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats.no_mask.png')
plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats_masked_adjusted.png',
                        trimming_mode='adjusted')
# Not implemented
# plotter.flen_strat_plot(output_path='/home/kraemers/temp/mbias_stats_masked_minimal.png',
#                         trimming_mode='minimal')
