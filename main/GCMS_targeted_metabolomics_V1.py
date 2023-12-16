#################################### load packages
import os
import sys
sys.path.append(os.getcwd())
from utils import functions
import yaml
import pandas as pd
import numpy as np

#################################### configure parameters
current_dir = os.getcwd()
config_path = os.path.join(current_dir, 'utils', 'config.yaml')

with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# necessary parameters
file_in = config['file_in']
input_table_type = config['input_table_type']
n_header_lines = config['n_header_lines']
file_out = config['file_out']
save_data_in_excel = config['save_data_in_excel']
save_plot = config['save_plot']
show_plot = config['show_plot']
target_column = config['target_column']
removing_samples_pattern = config['removing_samples_regular_expression']
rename_sample_from = config['rename_samples_regular_expression']['from']
rename_sample_to = config['rename_samples_regular_expression']['to']
sample_info_sequences = config['sample_info_sequences']

# batch correction parameters with defaults
batch_correction_method = config['batch_correction_method']
avoid_negative_values = config['avoid_negative_values']
batch_group = config['batch_group']
batch_size = config['batch_size']
batch_correction_loess_fraction = config['batch_correction_loess_fraction']

# outlier truncation parameters with defaults
max_standard_deviation_outlier = config['max_standard_deviation_outlier']

# missing value imputation parameters with defaults
missing_value_as_zero = config['missing_value_as_zero']
knn_neighbours = config['knn_neighbours']

# calculate coefficient of variation and Sample to RPO ratio with defaults
max_CV = config['max_CV']
min_ratio = config['min_ratio']

# plotting parameters with defaults
pca_scale_factor = config['pca_scale_factor']
global_text_size = config['global_text_size']
subplot_size = config['subplot_size']
pca_remove_control = config['pca_remove_control']
pca_remove_pattern = config['pca_remove_pattern']
check_metabolite = config['check_metabolite']

#################################### prepare tables

# parse table
if input_table_type == ['Agilent_targeted_metabolomics']:
    df = functions.targeted_metabolomics_agilent_parser(file_path=file_in, n_header=n_header_lines,
                                                        needed_table_columns=target_column)

    # pivot table
    df_wide = df.pivot(index=['Name', 'Data File', 'Type', 'Level', 'Acq. Date-Time', 'Metabolite'],
                       columns='variable',
                       values='value').reset_index()

# remove technical samples,
if removing_samples_pattern is not None:
    df_wide = functions.remove_samples(dataframe=df_wide, pattern=removing_samples_pattern)

# rename sample names
if rename_sample_from is not None:
    df_wide = functions.rename_samples(dataframe=df_wide, from_name=rename_sample_from, to_name=rename_sample_to)

# extract needed columns
df_conc = pd.pivot_table(df_wide, index='Name', columns='Metabolite', values=target_column[0])
df_conc = df_conc.astype('float')
df_conc = df_conc.dropna(axis=1, how='all')
df_conc = functions.reorder_dataframe(dataframe=df_conc, by='Run_number', sampleinfo=sample_info_sequences)
dict_conc = {'S0_Original': df_conc}

#################################### step 1: batch correction

# split data into batches
if batch_group == 'Split':
    df_conc = functions.batch_split(dataframe=df_conc, batch_size=batch_size, sampleinfo=sample_info_sequences)

df_conc.index.get_level_values('Split').value_counts()

# batch correction
if batch_correction_method == 'Default':
    dict_conc['S1_Batch_corrected'] = functions.batch_correct_interquartile(dataframe=df_conc,
                                                                            group=batch_group,
                                                                            sampleinfo=sample_info_sequences)
elif batch_correction_method == 'loess':
    functions.batch_correct_loess(dataframe=df_conc, group=batch_group, frac=batch_correction_loess_fraction,
                                  sampleinfo=sample_info_sequences, knn_neighbor=knn_neighbours)

#################################### step 2: outliers truncation

dict_conc['S2_outlier_truncated'] = functions.outlier_truncate(dataframe=dict_conc['S1_Batch_corrected'],
                                                               n_std=max_standard_deviation_outlier, log10=True)

#################################### step 3: missing value imputation
if missing_value_as_zero:
    dict_conc['S3_missing_value_as_zero'] = dict_conc['S2_outlier_truncated'].fillna(value=0)
else:
    dict_conc['S3_missing_value_imputed'] = functions.knn_impute(
    dataframe=dict_conc['S2_outlier_truncated'], n_neighbors=knn_neighbours)


#################################### step 4: calculating CV and Ratio
dict_cv = functions.compute_CV(data_dict=dict_conc, sampleinfo=sample_info_sequences, group=None, ratio_to='RPO')

if max_CV & min_ratio is not None:
    last_key = list(dict_cv.keys())[-1]
    df_last_pruned = dict_cv[last_key].loc[(dict_cv[last_key]['Sample'] <= max_CV) &
                                                       (dict_cv[last_key]['Ratio'] >= min_ratio), :]
    metabolites = df_last_pruned.index.values
    for key,df in dict_cv.items():
        dict_cv[key] = None
        dict_cv[key] = df.loc[metabolites,:]
    for key,df in dict_conc.items():
        dict_conc[key] = None
        dict_conc[key] = df.loc[:,metabolites]

#################################### step 5: save plots

functions.pca_plot(data_dict=dict_conc, group='Sample_type', sampleinfo=sample_info_sequences,
                   knn_neighbours=knn_neighbours, pca_scale_factor=pca_scale_factor,
                   global_text_size=global_text_size, subplot_size=subplot_size,
                   remove_pattern=pca_remove_pattern, remove_controls=False, save_plot=save_plot,
                   file_in=file_in, show_plot=show_plot)

functions.pca_plot(data_dict=dict_conc, group=batch_group, sampleinfo=sample_info_sequences,
                   knn_neighbours=knn_neighbours, pca_scale_factor=pca_scale_factor,
                   global_text_size=global_text_size, subplot_size=subplot_size,
                   remove_pattern=pca_remove_pattern, remove_controls=True, save_plot=save_plot,
                   file_in=file_in, show_plot=show_plot)

functions.scatter_plot(data_dict=dict_conc, file_in=file_in, sampleinfo=sample_info_sequences, row_height=2.5,
                       width2height_ratio=1.2, log10=False, group=batch_group, check_metabolite=check_metabolite,
                       alpha=0.5, size=1, y_value=target_column[0], show_plot=show_plot, save_plot=save_plot)

functions.CV_plot(data_dict=dict_cv, title_font_size=4, table_width=1.18, cell_font_size=2.5, save_plot=save_plot,
                  file_in=file_in, table_title_position=1, show_plot=show_plot)

#################################### step 6: save data to csv

if save_data_in_excel:
    for key,df in dict_conc.items():
        #key = 'S1_Batch_corrected'
        sample_info = functions.get_sample_info(df, columns=sample_info_sequences)
        index = list(df.index.names)
        index_name = [i for i in index if i not in sample_info_sequences + ['Split']]
        df.reset_index(inplace=True)
        df = sample_info.merge(df, on=index_name)
        df.set_index(index_name, inplace=True)
        dict_conc[key] = df

    items = list(dict_conc.items())

    # in case there is
    file_out = os.path.basename(file_out)
    dirname = os.path.dirname(file_in)
    dirname = os.path.join(dirname, 'output')
    file_out = os.path.join(dirname, file_out)
    if os.path.exists(file_out):
        raise KeyError('The output exists already!')

    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with pd.ExcelWriter(file_out, mode='w') as writer:
        sheetname = items[0][0]
        df = items[0][1]
        df.to_excel(writer, sheet_name=sheetname, index=True)

    for key, value in items[1:]:
        with pd.ExcelWriter(file_out, mode='a') as writer:
            value.to_excel(writer, sheet_name=key, index=True)

print('\nanalysis finished')