import pandas as pd

def humann_pathway(
    pathway_table: pd.DataFrame,
    strip_units_from_sample_ids: bool = True,
    destratify: bool = False) -> (pd.DataFrame, pd.DataFrame):

    pathway_table = pathway_table.reset_index()
    if destratify:
        pathway_table['unstratified'] = pathway_table.apply(
            lambda x: '|' not in x['feature-id'], axis=1)
        pathway_table = pathway_table[pathway_table['unstratified']]
        pathway_table = pathway_table.drop('unstratified', axis=1)

    # Generate the taxonomy result
    taxonomy = pathway_table[['feature-id']].copy()
    taxonomy['Taxon'] = taxonomy.apply(
        lambda x: x['feature-id'].replace('|', '; ').replace('.', '; '), axis=1)

    taxonomy = taxonomy.set_index('feature-id')
    taxonomy.index.name = 'Feature ID'


    # Generate the table
    pathway_table = pathway_table.set_index('feature-id')
    pathway_table = pathway_table.T
    pathway_table.index.name = 'sample-id'

    if strip_units_from_sample_ids:
        pathway_table = pathway_table.reset_index()

        pathway_table['sample-id'] = pathway_table.apply(
            lambda x: x['sample-id'].rsplit('_', 1)[0], axis=1)

        pathway_table = pathway_table.set_index('sample-id')


    return pathway_table, taxonomy
