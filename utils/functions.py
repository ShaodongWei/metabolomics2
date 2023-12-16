# needed packages: pandas, numpy, openpyxl, PyYAML, dateutil, matplotlib, sklearn, seaborn, random, os, statsmodels,


def convert_value(Series):
    """
    Convert the value to its potential type.
    """

    from dateutil.parser import parse
    try:
        # Try integer
        return Series.astype('int')
    except ValueError:
        pass

    try:
        # Try float
        return Series.astype('float')
    except ValueError:
        pass

    try:
        # Try date using dateutil.parser
        return Series.apply(lambda x: parse(x).date())  # returning just the date part
    except (ValueError, TypeError):
        pass

    return Series  # If all else fails, return the original string


def reorder_dataframe(dataframe, by, sampleinfo):
    # dataframe=df_conc
    # by='Run_number'
    # sampleinfo=sample_info_sequences

    '''
    reorder dataframe based on the information from sampleinfo (sample naming)

    '''
    import pandas as pd
    df = dataframe.copy()
    sample_info = get_sample_info(dataframe=df, columns=sampleinfo)
    sample_info[by] = convert_value(sample_info[by])
    sample_info = sample_info.sort_values(by=by)
    index = list(df.index.names)
    df.reset_index(inplace=True)
    merge_on = [i for i in df.columns if i in sample_info.columns]
    df = df.merge(sample_info[merge_on], on=merge_on, how='right')
    df.set_index(index, inplace=True)
    return df


def targeted_metabolomics_agilent_parser(file_path, n_header, needed_table_columns):
    '''
    This functions passes raw exce or csv table from Agilent & targeted & metabolomics
    '''
    import pandas as pd
    import numpy as np

    # input multiple-level header
    def get_file_extension(file):
        extension = file.split('.')[-1]
        return extension

    if get_file_extension(file=file_path) == 'csv':
        df = pd.read_csv(file_path, header=list(range(0, n_header)))
    elif get_file_extension(file=file_path) == 'xlsx':
        df = pd.read_excel(file_path, header=list(range(0, n_header)))

    # remove fist 2 columns
    if (df.iloc[:, 1] == "!").mean() >= 0.8:
        df = df.iloc[:, 2:]
    header = pd.Series(df.columns.get_level_values(0))  # make a new header

    # rename header
    metabolite = ""
    for i in range(len(header)):
        if 'Result' in header[i]:
            metabolite = header[i]
        elif 'Unname' in header[i]:
            header[i] = metabolite

    header = pd.Series([i.replace(" Results", "") for i in header])
    # header.value_counts()

    new_col = df.columns.get_level_values(1) + "_" + header

    # assign column with new names
    df.columns = new_col

    # melt data frame
    df_melt = pd.melt(df, id_vars=df.columns[0:5], value_vars=df.columns[5:])
    df_melt['Metabolite'] = df_melt['variable'].str.split("_").str[1].str.strip()
    df_melt['variable'] = df_melt['variable'].str.split('_').str[0].str.strip()
    df_melt.columns = df_melt.columns.str.replace('_', '')

    # remove duplicated rows
    check = df_melt['Name'] + df_melt['Metabolite']
    pos = check.duplicated()
    df_met = df_melt.loc[~pos, :]

    # check columns if having necessary measrements
    for measure in needed_table_columns:
        if measure not in needed_table_columns:
            raise KeyError(f"Your table does not have '{measure}' ")
        else:
            return df_melt


def remove_samples(dataframe, pattern, invert_results=False):
    '''
    This function removes (or keeps) samples based on a pattern list.
    '''

    import pandas as pd
    import numpy as np

    # check if patter is a list
    if not isinstance(pattern, list):
        pattern = [pattern]

    df = dataframe.copy()
    index = df.index.names

    # put df.index to columns if df have index
    if np.sum(pd.isna(list(index))) == 0:
        df.reset_index(inplace=True)

    mask = np.ones(len(df), dtype=bool)
    for i in pattern:
        pos = ~df['Name'].str.contains(i, regex=True, case=False)
        mask = mask & pos

    if not invert_results:
        df = df.loc[mask, :]
    else:
        df = df.loc[~mask, :]

    # put index back
    if np.sum(pd.isna(list(index))) == 0:
        df.set_index(index, inplace=True)

    return df


def rename_samples(dataframe, from_name, to_name):
    '''
    This function renames samples based on regular expression patterns.
    '''
    if not isinstance(from_name, list):
        raise KeyError('"from_name" has to be a list')
    if not isinstance(to_name, list):
        raise KeyError('"to_name" has to be a list')
    for i in range(len(from_name)):
        from_name_i = from_name[i]
        to_name_i = to_name[i]

        for col in dataframe.columns:
            try:
                if dataframe[col].str.contains(from_name_i).sum() > 0:
                    dataframe[col] = dataframe[col].str.replace(from_name_i, to_name_i, regex=True)
            except AttributeError:
                pass
    return dataframe


def get_sample_info(dataframe, columns):
    '''
    This function gets information from the sample name, delimited by '_'. 'columns' is the sequence of information.
    '''
    import pandas as pd
    import numpy as np
    if 'Name' in dataframe.index.names:
        tmp = pd.DataFrame({'Name': dataframe.index.get_level_values('Name')})
        tmp[columns] = tmp['Name'].str.split('_', expand=True)
    elif 'Name' in dataframe.columns:
        tmp = dataframe[['Name']].copy()
        tmp[columns] = tmp['Name'].str.split('_', expand=True)
    conditions = [
        tmp['Steno_label'].str.contains('^\d+$'),
        tmp['Steno_label'].str.contains('^cp', case=False),
        tmp['Steno_label'].str.contains('^rp', case=False)
    ]
    choices = ['Sample', 'CP', 'RPO']
    tmp['Sample_type'] = np.select(conditions, choices, default='Unknown')
    return tmp


def batch_split(dataframe, batch_size, sampleinfo, order_by='Run_number'):
    '''
    This function split samples into batches,  sized with 'batch_size', based on the order of 'order_by'.
    '''
    import pandas as pd
    import numpy as np

    sample_info = get_sample_info(dataframe=dataframe, columns=sampleinfo)

    for x, y in zip(sample_info['Name'], dataframe.index.get_level_values('Name')):
        if x != y:
            raise KeyError('sample_info Name different from dataframe Name')

    if order_by not in sampleinfo:
        raise KeyError('The order_by is not in sampleinfo')
    else:
        if 'Name' not in dataframe.index.names and 'Name' not in dataframe.columns:
            raise KeyError('The dataframe has no Name (index or column)')
        else:
            if 'Name' in dataframe.columns:
                dataframe.set_index('Name', inplace=True)
            else:
                dataframe_corrected = pd.DataFrame()
                batch_number = 1
                for batch in sample_info['Batch'].unique():
                    # batch = '20230526'

                    dataframe_batch = dataframe.loc[list(sample_info['Batch'] == batch), :]
                    if 'Split' in dataframe_batch.index.names:
                        dataframe_batch.reset_index(drop=True, level='Split', inplace=True)

                    sample_info_batch = sample_info.loc[sample_info['Batch'] == batch, :]
                    sample_info_batch.loc[:, order_by] = sample_info_batch[order_by].astype('int')
                    order = sample_info_batch.sort_values(order_by)['Name']
                    dataframe_batch_ordered = dataframe_batch.loc[order]
                    dataframe_batch_ordered['Run_order'] = range(0, dataframe_batch_ordered.shape[0])

                    split_number = int(np.round(dataframe_batch.shape[0] / int(batch_size)))
                    import math
                    number_in_split = math.ceil(dataframe_batch.shape[0] / split_number)
                    breaks = [i * number_in_split for i in range(0, split_number + 1)]
                    # breaks = [i for i in breaks if i != 0]

                    dataframe_split = pd.DataFrame(
                        pd.cut(dataframe_batch_ordered['Run_order'], bins=breaks, include_lowest=False, right=True))
                    dataframe_split = dataframe_split.astype('str')
                    dataframe_split.iloc[0, 0] = dataframe_split.iloc[1, 0]
                    dataframe_split.columns = ['Split']
                    # dataframe_split['Split'].value_counts() #this is same as Tommi's R code
                    dataframe_split['Split'] = ['Batch' + str(batch_number) + f'{i}' for i in dataframe_split['Split']]
                    dataframe_split.reset_index(inplace=True)
                    dataframe_batch_ordered.reset_index(inplace=True)
                    merge_on = [i for i in dataframe_batch_ordered.columns if i in dataframe_split.columns]
                    dataframe_batch_ordered = dataframe_batch_ordered.merge(dataframe_split, on=merge_on, how='left')
                    dataframe_batch_ordered.drop('Run_order', axis=1, inplace=True)
                    dataframe_batch_ordered.set_index(['Name', 'Split'], inplace=True)
                    batch_number = batch_number + 1
                    dataframe_corrected = pd.concat([dataframe_corrected, dataframe_batch_ordered], axis=0)
                return (dataframe_corrected)


def fillna_with_batch_statistic(dataframe, group, sampleinfo, method='Median'):
    '''
    This function replaces the NaN value in a dataframe with the statistic 'Median' or 'Mean'. The statistic is calculated
    with the 'group', which can be from the dataframe index or sampleinfo (sample_info_sequences).
    '''

    import warnings
    import pandas as pd
    import numpy as np
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    # Create a copy of the original DataFrame
    dataframe_copy = dataframe.copy()
    index_df = list(dataframe_copy.index.names)
    # dataframe_copy.reset_index(inplace=True)
    if group in dataframe_copy.index.names:
        dataframe_copy.reset_index(inplace=True)
    else:
        sample_info = get_sample_info(dataframe=dataframe_copy, columns=sampleinfo)
        sample_info = sample_info[['Name', group]]
        dataframe_copy.reset_index(inplace=True)
        merge_on = [i for i in dataframe_copy.columns if i in sample_info.columns]
        dataframe_copy = pd.merge(dataframe_copy, sample_info, on=merge_on)

    batches = dataframe_copy[group].unique()

    if method not in ['Median', 'Mean']:
        raise ValueError("method should be 'Median' or 'Mean'")

    for col in [i for i in dataframe_copy.columns if i not in ([group] + index_df)]:
        for bat in batches:
            subset = dataframe_copy.loc[dataframe_copy[group] == bat, col]
            statistic = subset.median() if method == 'Median' else subset.mean()
            if np.isnan(statistic):
                statistic = 0  # using zero, if median or mean is NA
            mask = (dataframe_copy[group] == bat) & (dataframe_copy[col].isna())
            dataframe_copy.loc[mask, col] = dataframe_copy.loc[mask, col].fillna(statistic)

    # Reset the index back to the original index
    dataframe_copy.set_index(list(set(index_df + [group])), inplace=True)

    return dataframe_copy

def batch_correct_interquartile(dataframe, group, sampleinfo, log10=False):
    '''
    This function gives same result as 'batch_correct_interquartile2', but coding differently. This methodology
    is based on Tommi's logic of coding. 'group' can be from the dataframe index or sampleinfo (sample_info_sequences)
    '''
    # do iqr correction first, then median correction, exactly same output compared to Tommi's method
    import pandas as pd
    import numpy as np

    # replace zero with NA
    dataframe[dataframe == 0] = np.nan

    # # replacing missing value
    # dataframe = fillna_with_batch_statistic(dataframe=dataframe, group=group, sampleinfo=sampleinfo, method='Median')

    df_batch_corrected = pd.DataFrame()

    def my_iqr(series):
        tmp_median_75 = np.nanpercentile(series.astype('float'), 75)
        tmp_median_25 = np.nanpercentile(series.astype('float'), 25)
        tmp_iqr = tmp_median_75 - tmp_median_25
        return tmp_iqr

    for i in dataframe.columns:
        # i='Valine'
        tmp_col = dataframe[[i]]

        if group in dataframe.index.names:
            # first correct based on IQR, but Tommis actually used mean of batch_IQRs
            global_iqr = np.mean(tmp_col.groupby(group).apply(my_iqr))
            global_median = np.mean(tmp_col.groupby(group).apply(np.nanmedian))

            if global_iqr == 0:
                raise KeyError('You only have one element in each group')
            else:
                tmp_col_iqr_corrected = pd.DataFrame()
                for j in tmp_col.index.get_level_values(group).unique():
                    # j='0-30'
                    tmp_iqr = tmp_col.xs(j, level=group, axis=0).copy()
                    batch_median = np.nanmedian(np.array(tmp_iqr).astype('float'))
                    batch_iqr = np.nanpercentile(np.array(tmp_iqr).astype('float'), 75) - np.nanpercentile(
                        np.array(tmp_iqr).astype('float'), 25)

                    if pd.isna(batch_iqr):
                        tmp_iqr_corrected_final = tmp_iqr.copy()  # if iqr is zero, we do nothing
                    else:
                        tmp_iqr_corrected_final = pd.DataFrame(
                            tmp_iqr.iloc[:, 0].astype('float') / batch_iqr * global_iqr).copy()
                    tmp_iqr_corrected_final.columns = [i]
                    tmp_iqr_corrected_final.reset_index(inplace=True)
                    tmp_iqr_corrected_final[group] = j
                    index = [col for col in tmp_iqr_corrected_final.columns if col != i]
                    tmp_iqr_corrected_final.set_index(index, inplace=True)
                    tmp_col_iqr_corrected = pd.concat([tmp_col_iqr_corrected, tmp_iqr_corrected_final], axis=0)

        elif group in sampleinfo:
            sample_info = get_sample_info(dataframe=tmp_col, columns=sampleinfo)
            tmp_col = tmp_col.copy()
            index = list(tmp_col.index.names)
            tmp_col = tmp_col.reset_index()
            merge_on = [i for i in tmp_col.columns if i in sample_info.columns]
            tmp_col = tmp_col.merge(sample_info[[group] + merge_on], on=merge_on)
            tmp_col.set_index(index + [group], inplace=True)

            # first correct based on IQR, but Tommis actually used mean of batch_IQRs
            global_iqr = np.mean(tmp_col.groupby(group).apply(my_iqr))
            global_median = np.mean(tmp_col.groupby(group).apply(np.nanmedian))
            if global_iqr == 0:
                raise KeyError('You only have one element in each group')
            else:
                tmp_col_iqr_corrected = pd.DataFrame()
                for j in sample_info[group].unique():
                    # j='20230523'
                    tmp_iqr = tmp_col.loc[sample_info[group].values == j, :].copy()
                    batch_median = np.nanmedian(np.array(tmp_iqr).astype('float'))
                    batch_iqr = np.nanpercentile(np.array(tmp_iqr).astype('float'), 75) - np.nanpercentile(
                        np.array(tmp_iqr).astype('float'), 25)
                    if pd.isna(batch_iqr) or batch_iqr == 0:
                        tmp_iqr_corrected_final = tmp_iqr.copy()  # if iqr is NaN, None or zero, we do nothing
                    else:
                        tmp_iqr_corrected_final = pd.DataFrame(
                            tmp_iqr.iloc[:, 0].astype('float') / batch_iqr * global_iqr).copy()
                    tmp_iqr_corrected_final.columns = [i]
                    tmp_iqr_corrected_final.reset_index(inplace=True)
                    tmp_iqr_corrected_final[group] = j
                    index = [col for col in tmp_iqr_corrected_final.columns if col != i]
                    tmp_iqr_corrected_final.set_index(index, inplace=True)
                    tmp_col_iqr_corrected = pd.concat([tmp_col_iqr_corrected, tmp_iqr_corrected_final], axis=0)

        else:
            raise KeyError('Your group is neither in sample_info_sequences nor in dataframe index')

        # next correct based on median
        tmp_col_median_corrected = pd.DataFrame()
        for k in tmp_col_iqr_corrected.index.get_level_values(group).unique():

            tmp_col_median = tmp_col_iqr_corrected.xs(k, level=group, axis=0).copy()
            batch_median = np.nanmedian(tmp_col_median)  # batch_median is actually batch_median*(global_iqr/batch_iqr)

            if not pd.isna(batch_median):
                tmp_median = tmp_col_median - batch_median + global_median
            tmp_median.columns = [i]
            tmp_median.reset_index(inplace=True)
            tmp_median[group] = k
            index = [col for col in tmp_median.columns if col != i]
            tmp_median.set_index(index, inplace=True)
            tmp_col_median_corrected = pd.concat([tmp_col_median_corrected, tmp_median], axis=0)

        df_batch_corrected = pd.concat([df_batch_corrected, tmp_col_median_corrected], axis=1)

    return df_batch_corrected


def batch_correct_interquartile2(dataframe, group, sampleinfo):
    '''
    dataframe = df_conc.copy()
    group = batch_group
    sampleinfo = sample_info_sequences
    '''

    '''
    This function gives same result as 'batch_correct_interquartile' (Tommi's method), but coding differently. This methodology
    is more intuitive in my opinion. 'group' can be from the dataframe index or sampleinfo (sample_info_sequences)
    '''
    # do iqr and median correction together. global median is real median, not mean(global group median as Tommi did), global
    # iqr is real iqr, not mean of global group iqr as Tommi did.

    # (x-median(batch)/iqr(batch), then x*iqr(global) + median(global)) does not work if median is zero !!!
    # the group has to be in the index, or in the sample_info_sequences

    import pandas as pd
    import numpy as np

    # replacing zero with missing value NaN
    dataframe[dataframe == 0] = np.nan

    # # replacing missing value
    # dataframe = fillna_with_batch_statistic(dataframe=dataframe, group=group, sampleinfo=sampleinfo, method='Median')

    df_batch_corrected = pd.DataFrame()
    for i in dataframe.columns:
        tmp_col = dataframe[[i]]
        global_median = np.nanmedian(np.array(tmp_col).astype('float'))
        global_median_75 = np.nanpercentile(np.array(tmp_col).astype('float'), 75)
        global_median_25 = np.nanpercentile(np.array(tmp_col).astype('float'), 25)
        global_iqr = global_median_75 - global_median_25

        if group in dataframe.index.names:
            tmp_col_corrected = pd.DataFrame()
            for j in tmp_col.index.get_level_values(group).unique():
                tmp_batch = tmp_col.xs(j, level=group, axis=0).copy()
                batch_median = np.nanmedian(np.array(tmp_batch).astype('float'))
                batch_iqr = np.nanpercentile(np.array(tmp_batch).astype('float'), 75) - np.nanpercentile(
                    np.array(tmp_batch).astype('float'), 25)

                if pd.isna(batch_iqr):
                    tmp_batch_corrected_final = tmp_batch.copy()  # if iqr is zero, we do nothing
                else:
                    tmp_batch_corrected1 = (tmp_batch.iloc[:, 0].astype('float') - batch_median) / batch_iqr
                    tmp_batch_corrected_final = pd.DataFrame(tmp_batch_corrected1 * global_iqr + global_median)
                tmp_batch_corrected_final.columns = [i]
                tmp_batch_corrected_final.reset_index(inplace=True)
                tmp_batch_corrected_final[group] = j
                index = [col for col in tmp_batch_corrected_final.columns if col != i]
                tmp_batch_corrected_final.set_index(index, inplace=True)
                tmp_col_corrected = pd.concat([tmp_col_corrected, tmp_batch_corrected_final], axis=0)
            df_batch_corrected = pd.concat([df_batch_corrected, tmp_col_corrected], axis=1)
        elif group in sampleinfo:
            sample_info = get_sample_info(dataframe=tmp_col, columns=sampleinfo)
            tmp_col_corrected = pd.DataFrame()
            for j in sample_info[group].unique():
                tmp_batch = tmp_col.loc[sample_info[group].values == j, :].copy()
                batch_median = np.nanmedian(np.array(tmp_batch).astype('float'))
                batch_iqr = np.nanpercentile(np.array(tmp_batch).astype('float'), 75) - np.nanpercentile(
                    np.array(tmp_batch).astype('float'), 25)
                if pd.isna(batch_iqr):
                    tmp_batch_corrected_final = tmp_batch.copy()  # if iqr is NaN, we do nothing
                else:
                    tmp_batch_corrected1 = (tmp_batch.iloc[:, 0].astype('float') - batch_median) / batch_iqr
                    tmp_batch_corrected_final = pd.DataFrame(tmp_batch_corrected1 * global_iqr + global_median)
                tmp_batch_corrected_final.columns = [i]
                tmp_batch_corrected_final.reset_index(inplace=True)
                tmp_batch_corrected_final[group] = j
                index = [col for col in tmp_batch_corrected_final.columns if col != i]
                tmp_batch_corrected_final.set_index(index, inplace=True)
                tmp_col_corrected = pd.concat([tmp_col_corrected, tmp_batch_corrected_final], axis=0)
            df_batch_corrected = pd.concat([df_batch_corrected, tmp_col_corrected], axis=1)
        else:
            df_batch_corrected = 'Your group is neither in sample_info_sequences nor in dataframe index'

    return df_batch_corrected

def batch_correct_interquartile3(dataframe, group, sampleinfo, avoid_negative_values):
    '''
    dataframe = df_conc.copy()
    group = batch_group
    sampleinfo = sample_info_sequences
    avoid_negative_values = True
    '''

    '''
    This function gives same result as 'batch_correct_interquartile2', but coding differently. This methodology
    is based on Tommi's logic of coding. 'group' can be from the dataframe index or sampleinfo (sample_info_sequences)
    '''
    # do iqr correction first, then median correction, exactly same output compared to Tommi's method
    import pandas as pd
    import numpy as np

    # replacing zero with missing value NaN
    #dataframe.replace({np.nan: 0}, inplace=True)
    dataframe[dataframe == 0] = np.nan

    # to avoid negative values after batch correction, we first do log10 transformation, then exponentiate it
    if avoid_negative_values:
        # replace NA with zero in order to do log10 transformation
        #dataframe.fillna(0, inplace=True)
        dataframe = dataframe.map(lambda x: np.log10(x+1))

    df_batch_corrected = pd.DataFrame()

    def my_iqr(series):
        tmp_median_75 = np.nanpercentile(series.astype('float'), 75)
        tmp_median_25 = np.nanpercentile(series.astype('float'), 25)
        tmp_iqr = tmp_median_75 - tmp_median_25
        return tmp_iqr

    for i in dataframe.columns:
        # i='Valine'
        tmp_col = dataframe[[i]]

        if group in dataframe.index.names:
            # first correct based on IQR, but Tommi actually used mean of batch_IQRs
            global_iqr = np.mean(tmp_col.groupby(group).apply(my_iqr))
            global_median = np.mean(tmp_col.groupby(group).apply(np.nanmedian))

            if global_iqr == 0:
                raise KeyError('You only have one element in each group')
            else:
                tmp_col_iqr_corrected = pd.DataFrame()
                for j in tmp_col.index.get_level_values(group).unique():
                    # j='0-30'
                    tmp_iqr = tmp_col.xs(j, level=group, axis=0).copy()
                    batch_median = np.nanmedian(np.array(tmp_iqr).astype('float'))
                    batch_iqr = np.nanpercentile(np.array(tmp_iqr).astype('float'), 75) - np.nanpercentile(
                        np.array(tmp_iqr).astype('float'), 25)

                    if pd.isna(batch_iqr):
                        tmp_iqr_corrected_final = tmp_iqr.copy()  # if iqr is zero, we do nothing
                    else:
                        tmp_iqr_corrected_final = pd.DataFrame(
                            tmp_iqr.iloc[:, 0].astype('float') / batch_iqr * global_iqr).copy()
                    tmp_iqr_corrected_final.columns = [i]
                    tmp_iqr_corrected_final.reset_index(inplace=True)
                    tmp_iqr_corrected_final[group] = j
                    index = [col for col in tmp_iqr_corrected_final.columns if col != i]
                    tmp_iqr_corrected_final.set_index(index, inplace=True)
                    tmp_col_iqr_corrected = pd.concat([tmp_col_iqr_corrected, tmp_iqr_corrected_final], axis=0)

        elif group in sampleinfo:
            sample_info = get_sample_info(dataframe=tmp_col, columns=sampleinfo)
            tmp_col = tmp_col.copy()
            index = list(tmp_col.index.names)
            tmp_col = tmp_col.reset_index()
            merge_on = [i for i in tmp_col.columns if i in sample_info.columns]
            tmp_col = tmp_col.merge(sample_info[[group] + merge_on], on=merge_on)
            tmp_col.set_index(index + [group], inplace=True)

            # first correct based on IQR, but Tommi actually used mean of batch_IQRs
            global_iqr = np.mean(tmp_col.groupby(group).apply(my_iqr))
            global_median = np.mean(tmp_col.groupby(group).apply(np.nanmedian))
            if global_iqr == 0:
                raise KeyError('You only have one element in each group')
            else:
                tmp_col_iqr_corrected = pd.DataFrame()
                for j in sample_info[group].unique():
                    # j='20230523'
                    tmp_iqr = tmp_col.loc[sample_info[group].values == j, :].copy()
                    batch_median = np.nanmedian(np.array(tmp_iqr).astype('float'))
                    batch_iqr = np.nanpercentile(np.array(tmp_iqr).astype('float'), 75) - np.nanpercentile(
                        np.array(tmp_iqr).astype('float'), 25)
                    if pd.isna(batch_iqr) or batch_iqr == 0:
                        tmp_iqr_corrected_final = tmp_iqr.copy()  # if iqr is NaN, None or zero, we do nothing
                    else:
                        tmp_iqr_corrected_final = pd.DataFrame(
                            tmp_iqr.iloc[:, 0].astype('float') / batch_iqr * global_iqr).copy()
                    tmp_iqr_corrected_final.columns = [i]
                    tmp_iqr_corrected_final.reset_index(inplace=True)
                    tmp_iqr_corrected_final[group] = j
                    index = [col for col in tmp_iqr_corrected_final.columns if col != i]
                    tmp_iqr_corrected_final.set_index(index, inplace=True)
                    tmp_col_iqr_corrected = pd.concat([tmp_col_iqr_corrected, tmp_iqr_corrected_final], axis=0)

        else:
            raise KeyError('Your group is neither in sample_info_sequences nor in dataframe index')

        # next correct based on median
        tmp_col_median_corrected = pd.DataFrame()
        for k in tmp_col_iqr_corrected.index.get_level_values(group).unique():

            tmp_col_median = tmp_col_iqr_corrected.xs(k, level=group, axis=0).copy()
            batch_median = np.nanmedian(tmp_col_median)  # batch_median is actually batch_median*(global_iqr/batch_iqr)

            if not pd.isna(batch_median):
                tmp_median = tmp_col_median - batch_median + global_median
            tmp_median.columns = [i]
            tmp_median.reset_index(inplace=True)
            tmp_median[group] = k
            index = [col for col in tmp_median.columns if col != i]
            tmp_median.set_index(index, inplace=True)
            tmp_col_median_corrected = pd.concat([tmp_col_median_corrected, tmp_median], axis=0)

        df_batch_corrected = pd.concat([df_batch_corrected, tmp_col_median_corrected], axis=1)


    if avoid_negative_values:
        df_batch_corrected.map(lambda x: 10**x - 1)
        # to convert zeros back to missing value
        #df_batch_corrected.replace({0: np.nan}, inplace=True)

    return df_batch_corrected

def batch_correct_loess(dataframe, group, frac, sampleinfo, knn_neighbor):
    '''
    This function does LOESS (Locally Estimated Scatterplot Smoothing) batch correction.
    'group' is 'x' or a variable that you want to smooth 'y' value along with.
    'group' can be from dataframe index or 'sampleinfo'.
    'frac' is a sensitive parameter (0 to 1) that controls the extent of smooth.
    'knn_neighbor' is the number of neighbor values to infer the missing value.
    '''

    # group is from either the sample_info_sequences or dataframe index
    import pandas as pd
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess

    # imputing missing value
    dataframe = knn_impute(dataframe=dataframe, n_neighbors=knn_neighbor)
    # dataframe = fillna_with_batch_statistic(dataframe=dataframe, group=group, sampleinfo=sampleinfo, method='Median')

    df_batch_corrected = pd.DataFrame()
    if group in dataframe.index.names:
        batch = dataframe.index.get_level_values(group)
        batch_test = batch[0]
        # Check if batch is already a numerical type (integer or float)
        if isinstance(batch_test, (int, float)):
            # If batch is a number, use it directly
            numerical_batch = batch
        elif isinstance(batch_test, str) and batch_test.replace('.', '', 1).isdigit():
            # If batch is a string representing an integer or float, convert it to float
            numerical_batch = batch.astype('float')
        else:
            batch = pd.Series(batch).astype('category').cat.codes

        for metabolite in dataframe.columns:
            # metabolite='Valine'
            y = dataframe[metabolite]
            # Fit the LOWESS model
            z = lowess(y, batch, frac=frac)

            # create a dictionary
            dict_y = {xvalue: yvalue for xvalue, yvalue in z}

            # Create batch_effect in the original order of the batch codes
            batch_effect = [dict_y[batch_code] for batch_code in batch]

            # Subtract the batch effect from the original data
            df_batch_corrected[metabolite] = y - batch_effect

        return df_batch_corrected
    elif group in sampleinfo:

        sample_info = get_sample_info(dataframe=dataframe, columns=sampleinfo)
        batch = sample_info[group]
        batch_test = batch[0]
        # Check if batch is already a numerical type (integer or float)
        if isinstance(batch_test, (int, float)):
            # If batch is a number, use it directly
            batch = batch
        elif isinstance(batch_test, str) and batch_test.replace('.', '', 1).isdigit():
            # If batch is a string representing an integer or float, convert it to float
            batch = batch.astype('float')
        else:
            batch = sample_info[group].astype('category').cat.codes

        df_batch_corrected = pd.DataFrame()
        for metabolite in dataframe.columns:
            # metabolite='Valine'
            y = dataframe[metabolite]
            # Fit the LOWESS model
            z = lowess(y, batch, frac=frac)

            # create a dictionary
            dict_y = {xvalue: yvalue for xvalue, yvalue in z}

            # Create batch_effect in the original order of the batch codes
            batch_effect = [dict_y[batch_code] for batch_code in batch]

            # Subtract the batch effect from the original data
            df_batch_corrected[metabolite] = y - batch_effect
        sample_info = get_sample_info(dataframe=df_batch_corrected, columns=sampleinfo)
        index = list(df_batch_corrected.index.names)
        df_batch_corrected.reset_index(inplace=True)
        merge_on = [i for i in df_batch_corrected.columns if i in sample_info.columns]
        df_batch_corrected = df_batch_corrected.merge(sample_info[[group] + merge_on], on=merge_on)
        df_batch_corrected.set_index(index + [group], inplace=True)
        return df_batch_corrected

    import os


def pca_plot(data_dict, group, sampleinfo, knn_neighbours, file_in, pca_scale_factor=1,
             global_text_size=8, subplot_size=5, remove_controls=False, remove_pattern=None,
             save_plot=True, show_plot=False):
    '''
    This function plots PCA for data collected in each step.
    'group' is used for color and ellipse, which can be from dataframe index or 'sampleinfo'.
    'knn_neighbours' is the number of neighbor values to infer the missing value.
    'file_out' is the saved output excel file, saved pca plots will stay with the excel file.
    'pca_scale_factor' is scaling factor between 0 and 1 that scales the size of the ellipse.
    'global_text_size' controls the overall text size.
    'subplot_size' is both the figure width and height.
    'remove_controls' controls if to remove controls based on 'remove_pattern'.
    'save_plot' controls if to save the pca plot as pdf.
    'show_plot' controls if to pops up plot on the windows.
    '''
    import pandas as pd
    import numpy as np
    import os
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from sklearn.decomposition import PCA

    # if or not remove CP, RPO
    if remove_controls and remove_pattern is not None:
        data_dict_updated = data_dict.copy()
        for key in data_dict.keys():
            tmp = data_dict[key].copy()
            tmp = remove_samples(dataframe=tmp, pattern=remove_pattern)
            data_dict_updated[key] = tmp
    else:
        data_dict_updated = data_dict.copy()

    # check where the group is
    presence_index = []
    for i in data_dict_updated.keys():
        # i = 'S1_Split'
        tmp = data_dict_updated[i]
        if group in tmp.index.names:
            presence_index.append(i)
    if len(presence_index) > 0:
        presence_group = 'index'

        # add group to index if one dataframe does not have e.g. S0_Original
        key_miss_group = []
        for i in data_dict_updated.keys():
            tmp = data_dict_updated[i]
            if group not in tmp.index.names:
                key_miss_group.append(i)

        key_not_miss_group = [i for i in data_dict_updated.keys() if i not in key_miss_group]
        for key in key_miss_group:
            # key = 'S1_Original'
            tmp = data_dict_updated[key].copy()
            other = key_not_miss_group[0]
            data_dict_updated_other = data_dict_updated[other]
            index_df = pd.DataFrame(data_dict_updated_other.index.tolist(), columns=data_dict_updated_other.index.names)
            tmp = data_dict_updated[key].copy()
            tmp.reset_index(inplace=True)
            merge_on = [i for i in tmp.columns if i in index_df.columns]
            tmp = tmp.merge(index_df, on=merge_on, how='left')
            tmp.set_index(list(index_df.columns), inplace=True)
            data_dict_updated[key] = tmp
    else:
        presence_group = 'sampleinfo'

        # Determine the number of plots
    num_plots = len(data_dict_updated)

    # Decide on the layout of the subplots based on the number of plots
    if num_plots > 2:
        num_cols = 2
        num_rows = int(np.ceil(num_plots / num_cols))
    else:
        num_cols = num_plots
        num_rows = 1

    fig_width = num_cols * subplot_size
    fig_height = num_rows * subplot_size

    if not show_plot:
        plt.ioff()  # Turn off interactive mode

    # Create a grid of subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(fig_width, fig_height), squeeze=False)
    axes = axes.flatten()

    # Loop through each key-value pair and each axis
    for ax, (key, value) in zip(axes, data_dict_updated.items()):

        if presence_group == 'index':
            group_info = pd.DataFrame(value.index.get_level_values(group))
        elif presence_group == 'sampleinfo':
            sample_info = get_sample_info(dataframe=value, columns=sampleinfo)
            group_info = sample_info[[group]]

        # impute missing value
        df = knn_impute(dataframe=value, n_neighbors=knn_neighbours, log10=True)

        # Perform PCA
        pca_out = PCA().fit_transform(df)
        pca_out = pd.DataFrame(pca_out[:, 0:2], columns=['PC1', 'PC2'])
        pca_out[group] = group_info

        # for plotting
        groups = pca_out[group].unique()
        colors = sns.color_palette("Set1", n_colors=len(groups))

        # Adjust legend
        
        plt.rcParams['font.size'] = global_text_size
        sns.scatterplot(data=pca_out, x='PC1', y='PC2', hue=group, palette=colors, ax=ax)
        ax.tick_params(axis='both', labelsize=global_text_size)
        ax.set_xlabel('PC1', fontsize=global_text_size + 2)
        ax.set_ylabel('PC2', fontsize=global_text_size + 2)
        ax.set_aspect('equal', 'box')
        ax.legend(title=group)

        for i, batch in enumerate(groups):
            group_data = pca_out[pca_out[group] == batch]
            mean = group_data[['PC1', 'PC2']].mean().values
            cov = np.cov(group_data['PC1'], group_data['PC2'])
            eigenvalues, eigenvectors = np.linalg.eig(cov)
            angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
            width, height = 2 * np.sqrt(eigenvalues) * pca_scale_factor
            ellipse = patches.Ellipse(mean, width, height, angle=angle, edgecolor=colors[i], facecolor='none')
            ax.add_patch(ellipse)
        ax.set_title(f'{key}')

    # Turn off the axis for any remaining subplots
    for i in range(num_plots, num_rows * num_cols):
        axes[i].axis('off')

    if save_plot:
        dirname = os.path.dirname(file_in)
        dirname = os.path.join(dirname, 'output')
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        filename = get_unique_filename(directory=dirname, base_name='pca_plot', extension='.pdf')
        plt.savefig(filename, format='pdf')
    if show_plot:
        plt.tight_layout()
        plt.show()

def scatter_plot(data_dict, file_in, sampleinfo, y_value, row_height=3, width2height_ratio=1.2, log10=True,
                 group='Run_number', check_metabolite='', alpha=0.5, size=1,
                 show_plot=False, save_plot=True, offset=0):
    '''
    This function plots a scatterplot with 5 metabolites as rows, each step of processing as columns.
    'group' is the 'x' axis in each subplot, axis will be ordered based on it.
    'check' is a list or None or '', used to check specific metabolites, e.g. ['a','b'].
    '''

    # if specify group, it is ordered and grouped based on group
    # if group is Split, ordered by the Run_number

    import pandas as pd
    import numpy as np
    import random
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os

    # add zeros with an offset
    for key,df in data_dict.items():
        df.replace(to_replace=offset, value=0, inplace=True)
        data_dict[key] = df

    if check_metabolite:
        first_key = list(data_dict.keys())[0]
        if check_metabolite not in data_dict[first_key].columns:
            raise KeyError('check_metabolite is not in the table !!!' )
        batch = [check_metabolite]
    else:
        # take turns to show 5 metabolite at a time
        first_key = list(data_dict.keys())[0]
        elements = list(data_dict[first_key].columns)
        random.shuffle(elements)
        # Divide elements into batches of 5
        batches = [elements[i:i + 5] for i in range(0, len(elements), 5)]

        # Define a generator function that iterates through batches
        def get_next_batch():
            while True:  # Infinite loop to keep restarting the iteration
                for batch in batches:
                    yield batch

        # Create a generator object
        next_batch_generator = get_next_batch()

        # Function to call
        def get_batch():
            return next(next_batch_generator)

        batch = get_batch()  # Returns the first batch

    df_ot = pd.DataFrame()
    for key in list(data_dict.keys())[::-1]:
        # key = 'S0_Original'
        if key != list(data_dict.keys())[0]:
            #key = 'S1_Batch_corrected'
            #key = 'S2_outlier_truncated'
            subset = data_dict[key][batch].copy()
            index_name = subset.index.names
            index_merge = pd.DataFrame(list(subset.index), columns=index_name)
            subset.reset_index(index_name, inplace=True)
            subset = pd.melt(subset, id_vars=index_name, var_name=['variable'], value_name='value').copy()
            subset['Data'] = key
            subset[sampleinfo] = subset['Name'].str.split('_', expand=True)
            df_ot = pd.concat([df_ot, subset], axis=0)
        else:
            #key = 'S0_Original'
            subset = data_dict[key][batch].copy()
            index = subset.index.names
            subset.reset_index(index, inplace=True)
            subset = pd.melt(subset, id_vars=index, var_name=['variable'], value_name='value').copy()
            subset['Data'] = key
            subset[sampleinfo] = subset['Name'].str.split('_', expand=True)
            subset = subset.merge(index_merge, on='Name', how='left')
            df_ot = pd.concat([df_ot, subset], axis=0)

    if log10:
        # Apply the logarithmic transformation to the 'value' column
        df_ot['log_value'] = df_ot['value'].apply(lambda x: 0 if x <= 0 else np.log(x))
    else:
        df_ot['log_value'] = df_ot['value']

    # re-order dataframe
    df_ot_sorted = df_ot.sort_values(by=group)

    # order columns
    col_order = list(df_ot_sorted['Data'].unique())

    def extract_number(element):
        return int(element.split('S')[1].split('_')[0])

    col_order = sorted(col_order, key=extract_number)

    # Creating the FacetGrid and plotting
    g = sns.FacetGrid(df_ot_sorted, row="variable", col="Data", col_order=col_order,
                      height=row_height, aspect=width2height_ratio, sharey=False)

    def convert_categorical_to_numeric(series):
        if series.dtype.name == 'category':
            return series.cat.codes
        elif series.dtype == 'object':
            # Assuming the column is string type
            return series.astype('category').cat.codes
        else:
            # Column is already numeric
            return series

    df_ot_sorted['numeric_group'] = convert_categorical_to_numeric(df_ot_sorted[group])

    # Mapping the scatter plot to the FacetGrid with log-transformed y-values
    g.map_dataframe(sns.scatterplot, x=group, y="log_value", alpha=alpha, size=size)
    g.map_dataframe(sns.regplot, x='numeric_group', y="log_value", scatter=False, lowess=True, color='red')

    # Adding y-label with actual variable names to the left
    for metabolite, ax in zip(g.row_names, g.axes[:, 0]):
        # print(metabolite, ax)
        ax.set_ylabel(ylabel=metabolite, fontsize=10)

    # Adding y-label with actual variable names to the right and removing tick labels
    for ax, title in zip(g.axes[:, -1], g.row_names):
        right_axis = ax.twinx()
        right_axis.spines['right'].set_visible(False)
        right_axis.spines['top'].set_visible(False)
        right_axis.set_ylabel(title, rotation=-90, va='center', fontsize=10)
        right_axis.set_yticks([])

    # Manually setting all subplot titles (removing the vertical line)
    for row_axes in g.axes:
        for ax, col_name in zip(row_axes, g.col_names):
            ax.set_title(col_name, fontsize=10)

    # disable the axis title from 2nd row
    for ax_row in g.axes[1:, :]:
        for ax in ax_row:
            ax.set_title('')

    for ax_row in g.axes:
        for ax in ax_row:
            ax.set_xticks([])

    # Setting the axis label size
    if log10:
        g.set_axis_labels(group, 'log10 of ' + y_value, fontsize=10)
    else:
        g.set_axis_labels(group, y_value, fontsize=10)

    # g.fig.suptitle('Scatter Plot with metabolites and correction methods', fontsize=12, y=0.98)
    g.tight_layout()

    for row_axes in g.axes:
        # Determine the y-limits for this row based on the data in the first subplot
        ylims = row_axes[0].get_ylim()
        # Set the y-limits for all subplots in this row to the determined limits
        for ax in row_axes:
            ax.set_ylim(ylims)

    # Displaying the plot
    if show_plot:
        plt.show()

    if save_plot:
        dirname = os.path.dirname(file_in)
        dirname = os.path.join(dirname, 'output')
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        filename = get_unique_filename(directory=dirname, base_name='scatter_plot', extension='.pdf')
        plt.savefig(filename, format='pdf')


def outlier_truncate(dataframe, n_std, log10=True):
    '''
    This function imputes outliers, any value > median + max_sd_allowed or < median - max_sd_allowed (outlier)
    will be assigned the max or min allowed value.

    '''

    import pandas as pd
    import numpy as np

    dataframe = dataframe.astype('float64')

    if log10:
        dataframe_updated = dataframe.apply(np.log10).copy()
    else:
        dataframe_updated = dataframe.copy()

    median = dataframe_updated.apply(np.nanmedian)
    sd = dataframe_updated.apply(np.nanstd, ddof=1)
    sd_allowed = sd * n_std

    max_outlier = median + sd_allowed
    min_outlier = median - sd_allowed

    def truncate(dataframe_series):
        series_name = dataframe_series.name
        max_allowed_value = max_outlier[series_name]
        min_allowed_value = min_outlier[series_name]

        dataframe_series[dataframe_series > max_allowed_value] = max_allowed_value
        dataframe_series[dataframe_series < min_allowed_value] = min_allowed_value

        return dataframe_series

    dataframe_updated = dataframe_updated.apply(truncate)
    if log10:
        dataframe_updated = 10 ** dataframe_updated
    return dataframe_updated


def knn_impute(dataframe, n_neighbors, log10=True):
    '''
    This function uses KNNImputer to impute missing value, it gives different slight results than R package impute::impute.knn
    '''

    # this function imputation follow steps: log10, (x - mean)/sd, knn, x*sd + mean, 10**x
    # same as what Tommi did
    import numpy as np
    import pandas as pd
    from sklearn.impute import KNNImputer

    if log10:
        df_updated = dataframe.apply(np.log10)

    df_updated[np.isinf(df_updated)] = np.nan  # replace inf with nan
    df_mean = df_updated[np.isfinite(df_updated)].mean(skipna=True)  # exlude inf and nan
    df_sd = df_updated[np.isfinite(df_updated)].std(skipna=True)  # exclude inf and nan
    df_scaled = (df_updated - df_mean) / df_sd

    knn = KNNImputer(n_neighbors=n_neighbors, weights='uniform', metric='nan_euclidean')
    df_knn = pd.DataFrame(knn.fit_transform(df_scaled))
    df_knn.columns = df_scaled.columns
    df_knn.index = df_scaled.index

    df_original = df_knn * df_sd + df_mean

    if log10:
        df_original = 10 ** df_original

    return df_original


def compute_CV(data_dict, sampleinfo, group=None, ratio_to='RPO'):
    '''
    This function calculates CV: (std/mean)*100. 'Ratio' is 'Sample/RPO'.
    '''

    import numpy as np
    import pandas as pd

    data_dict_processed = {}

    def cv(dataframe_col):
        return np.nanstd(dataframe_col, ddof=1) * 100 / np.nanmean(dataframe_col)

    for key, df in data_dict.items():
        sample_info = get_sample_info(dataframe=df, columns=sampleinfo)
        if group is not None:
            if group in sample_info.columns:
                group_series = list(sample_info[group])
            if group in df.index.names:
                group_series = list(df.index.get_level_values(group))
            grouped = df.groupby(group_series)
            df_cv = grouped.agg(cv)
            df_cv = df_cv.transpose()
            df_cv = df_cv.sort_values(by='Ratio', ascending=False)

        else:
            cp = remove_samples(dataframe=df, pattern='_CP_', invert_results=True)
            rpo = remove_samples(dataframe=df, pattern='_RPO_', invert_results=True)
            sample = remove_samples(dataframe=df, pattern=['_RPO_', '_CP_'])
            cp, rpo, sample = cp.apply(cv), rpo.apply(cv), sample.apply(cv)
            df_cv = pd.concat([sample, cp, rpo], axis=1)
            df_cv.columns = ['Sample', 'CP', 'RPO']
            df_cv['Ratio'] = df_cv['Sample'] / df_cv['RPO']
            df_cv = df_cv.sort_values(by='Ratio', ascending=False)

        data_dict_processed[key] = df_cv

    # order tables
    key_last = list(data_dict_processed.keys())[-1]
    metabolites_order = pd.DataFrame(data_dict_processed[key_last].index, columns=['Metabolite'])

    data_dict_ordered = {}
    for key, df in data_dict_processed.items():
        if key != key_last:
            df = metabolites_order.merge(df.reset_index(), left_on='Metabolite', right_on='index').drop(
                columns='index').set_index('Metabolite')
            data_dict_ordered[key] = df
        else:
            data_dict_ordered[key] = df

    return data_dict_ordered


def CV_plot(data_dict, file_in, table_width=1, cell_height=0.03, title_font_size=4, cell_font_size=3,
            table_title_position=1, save_plot=True, show_plot=False):
    '''
    This function plots the CV tables.
    'file_out' determines where to save the cv plot.
    'table_width' is the width for each sub table.
    'table_title_relative_height' controls the height of table title, it can be a negative value.
    'save_plot' determines if to save the plot.
    'show_plot' determines if to show the plot on screen.
    '''

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import os

    num_plot = len(data_dict)

    # Decide on the layout of the subplots
    if num_plot > 4:
        num_cols = 4
        num_rows = int(np.ceil(num_plot / num_cols))
    else:
        num_cols = num_plot
        num_rows = 1

    if not show_plot:
        plt.ioff()  # Turn off interactive mode

    fig, axes = plt.subplots(num_rows, num_cols)
    axes = axes.flatten()

    for ax, (key, value) in zip(axes, data_dict.items()):
        df = value.map(lambda x: format(x, '.2f'))

        # Plot table
        index_name = 'Metabolites'
        all_columns = [index_name] + df.columns.tolist()
        table = ax.table(cellText=df.reset_index().values, colLabels=all_columns, cellLoc='center', loc='center')

        ax.axis('off')
        ax.set_title(key, y=table_title_position, fontsize=title_font_size)

        # Adjust cell width based on text content and normalize based on table_width
        padding = 0
        col_widths = [max(len(col), max(map(len, df[col].astype(str).tolist()))) + padding for col in df.columns]
        index_width = max(map(len, df.index.astype(str).tolist())) + padding
        all_widths = [index_width] + col_widths

        # Normalize cell widths based on the table_width
        total_text_based_width = sum(all_widths)
        normalized_widths = [(width / total_text_based_width) * table_width for width in all_widths]

        # Set cell widths and heights
        for (i, j), cell in table.get_celld().items():
            cell.set_width(normalized_widths[j])
            cell.set_height(cell_height)

            # Colors for alternating rows
            colors_dict = {
                'lightergray': (0.9, 0.9, 0.9),
                'white': 'white'
            }

            cell.set_linewidth(0)

            # List of colors for ease of indexing
            colors_list = [colors_dict['lightergray'], colors_dict['white']]

            if i == 0:  # If it's the header row
                cell.get_text().set_weight('bold')
            else:
                cell.set_facecolor(colors_list[i % 2])

        table.auto_set_font_size(False)
        table.set_fontsize(cell_font_size)

    if save_plot:
        dirname = os.path.dirname(file_in)
        dirname = os.path.join(dirname, 'output')
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        file_name = get_unique_filename(directory=dirname, base_name='CV_plot', extension='.pdf')
        plt.savefig(file_name, format='pdf')
    if show_plot:
        plt.show()


def get_unique_filename(directory, base_name='pca_plot', extension=".pdf"):
    '''
    This function avoids overwritting existing .pdf files. When it detects a pdf with 'base_name', it adds 1 to the end of it.
    '''
    import os

    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter files by the specified extension
    files_with_extension = [f for f in all_files if f.endswith(extension)]

    # Find a unique name for the new file
    counter = 1
    filename = f"{base_name}{counter}{extension}"

    while filename in files_with_extension:
        counter += 1
        filename = f"{base_name}{counter}{extension}"

    return os.path.join(directory, filename)
