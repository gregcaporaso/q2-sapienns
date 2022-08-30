import pandas as pd


def _humann(table, strip_units_from_sample_ids, destratify):

    table = table.reset_index()
    if destratify:
        table['unstratified'] = table.apply(
            lambda x: '|' not in x['feature-id'], axis=1)
        table = table[table['unstratified']]
        table = table.drop('unstratified', axis=1)

    # Generate the taxonomy result
    taxonomy = table[['feature-id']].copy()
    taxonomy['Taxon'] = taxonomy.apply(
        lambda x: x['feature-id'].replace('|', '; ').replace('.', '; '),
        axis=1)

    taxonomy = taxonomy.set_index('feature-id')
    taxonomy.index.name = 'Feature ID'

    # Generate the table
    table = table.set_index('feature-id')
    table = table.T
    table.index.name = 'sample-id'

    if strip_units_from_sample_ids:
        table = table.reset_index()

        table['sample-id'] = table.apply(
            lambda x: x['sample-id'].rsplit('_', 1)[0], axis=1)

        table = table.set_index('sample-id')

    return table, taxonomy


def humann_pathway(
        pathway_table: pd.DataFrame,
        strip_units_from_sample_ids: bool = True,
        destratify: bool = False) -> (pd.DataFrame, pd.DataFrame):
    return _humann(pathway_table, strip_units_from_sample_ids, destratify)


def humann_genefamily(
        genefamily_table: pd.DataFrame,
        strip_units_from_sample_ids: bool = True,
        destratify: bool = False) -> (pd.DataFrame, pd.DataFrame):
    return _humann(genefamily_table, strip_units_from_sample_ids, destratify)
