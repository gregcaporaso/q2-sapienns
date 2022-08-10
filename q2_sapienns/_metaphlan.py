import pandas as pd

def metaphlan_taxon(
    stratified_table: pd.DataFrame, level: int)\
    -> (pd.DataFrame, pd.DataFrame):

    # Add a column indicating the number of levels contained in each
    # feature id.
    stratified_table['n levels'] = stratified_table.apply(
        lambda x: len(x['NCBI_tax_id'].split('|')), axis=1)

    # Drop features where number of levels is not equal to what was requested
    # by the user to "de-stratify" the table.
    table = stratified_table[stratified_table['n levels'] == level]
    if table.shape[0] == 0:
        raise ValueError('No features contained exactly %d taxonomic levels.' % level)

    # Generate the taxonomy result
    taxonomy = table['NCBI_tax_id']
    taxonomy = taxonomy.reset_index()
    taxonomy['Taxon'] = taxonomy.apply(
        lambda x: x['feature-id'].replace('|', '; '), axis=1)
    taxonomy = taxonomy.drop('feature-id', axis=1)
    taxonomy = taxonomy.set_index('NCBI_tax_id')
    taxonomy.index.name = 'Feature ID'

    table = table.reset_index()
    table = table.drop(['feature-id', 'n levels'], axis=1)
    table = table.set_index('NCBI_tax_id')
    table.index.name = 'sample-id'

    return table.T, taxonomy
