# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def metaphlan_taxon(
        stratified_table: pd.DataFrame, level: int)\
        -> (pd.DataFrame, pd.DataFrame):

    stratified_table = stratified_table.reset_index()

    # Add a column indicating the number of levels contained in each
    # feature id.
    stratified_table['n levels'] = stratified_table.apply(
        lambda x: len(x['feature-id'].split('|')), axis=1)

    # Drop features where number of levels is not equal to what was requested
    # by the user to "de-stratify" the table.
    table = stratified_table[stratified_table['n levels'] == level]
    if table.shape[0] == 0:
        raise ValueError('No features contained exactly %d taxonomic '
                         'levels.' % level)

    # Generate the taxonomy result
    taxonomy = table['feature-id'].to_frame()
    taxonomy['Taxon'] = taxonomy.apply(
        lambda x: x['feature-id'].replace('|', '; '), axis=1)
    taxonomy = taxonomy.set_index('feature-id')
    taxonomy.index.name = 'Feature ID'

    table = table.drop(['NCBI_tax_id', 'n levels'], errors='ignore', axis=1)
    table = table.set_index('feature-id')
    table = table.T
    table.index.name = 'sample-id'

    return table, taxonomy


def frequency(table: pd.DataFrame, target_freq: int = 100000) -> pd.DataFrame:
    target_freq /= 100
    result = table * target_freq
    result = result.round(0)
    result = result.astype(int)
    return result
