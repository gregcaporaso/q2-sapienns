from qiime2.plugin import (Plugin, SemanticType, TextFileFormat, model,
                           ValidationError, Citations, Int, Range)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy

import q2_sapienns

import pandas as pd

plugin = Plugin(
    name='sapienns',
    version=q2_sapienns.__version__,
    website='https://qiime2.org',
    user_support_text='https://forum.qiime2.org',
    package='q2_sapienns'
)

MergedMetaphlanTable = SemanticType('MergedMetaphlanTable')

plugin.register_semantic_types(MergedMetaphlanTable)

class MergedMetaphlanTableFormat(TextFileFormat):

    def _equal_number_of_columns(self, n_lines):
        with self.open() as fh:
            header_line = fh.readline()
            n_header_fields = len(header_line.split('\t'))
            if n_header_fields < 2:
                raise ValidationError(
                    'No sample columns appear to be present.')
            for idx, line in enumerate(fh, 2):
                if n_lines is not None and idx > n_lines + 1:
                    break
                n_fields = len(line.split('\t'))
                if n_fields != n_header_fields:
                    raise ValidationError(
                        'Number of columns on line %d is inconsistent with '
                        'the header line.' % line)

    def _validate_(self, level):
        level_to_n_lines = {'min': 5, 'max': None}
        self._equal_number_of_columns(level_to_n_lines[level])

MergedMetaphlanTableDirectoryFormat = model.SingleFileDirectoryFormat(
    'MergedMetaphlanTableDirectoryFormat', 'table.tsv', MergedMetaphlanTableFormat)

plugin.register_formats(MergedMetaphlanTableFormat,
                        MergedMetaphlanTableDirectoryFormat)

plugin.register_semantic_type_to_format(MergedMetaphlanTable,
                                        MergedMetaphlanTableDirectoryFormat)


@plugin.register_transformer
def _1(ff: MergedMetaphlanTableFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', header=0, index_col=0,
                         comment='#')
    result.index.name = 'feature-id'
    return result

citations = Citations.load('citations.bib', package='q2_sapienns')

def select_stratum(stratified_table: pd.DataFrame, level: int = None)\
    -> (pd.DataFrame, pd.DataFrame):

    if level is not None:
        # Add a column indicating the number of levels contained in each
        # feature id.
        stratified_table['n levels'] = stratified_table.apply(
            lambda x: len(x['NCBI_tax_id'].split('|')), axis=1
        )

        # Drop features where number of levels is not equal to what was requested
        # by the user.
        stratified_table = stratified_table[stratified_table['n levels'] == level]
        if stratified_table.shape[0] == 0:
            raise ValueError('No features contained exactly %d taxonomic levels.' % level)
        columns_to_drop = ['feature-id', 'n levels']
    else:
        columns_to_drop = ['feature-id']

    # Generate the taxonomy result
    taxonomy = stratified_table['NCBI_tax_id']
    taxonomy = taxonomy.reset_index()
    taxonomy['Taxon'] = taxonomy.apply(
        lambda x: x['feature-id'].replace('|', '; '), axis=1
    )
    taxonomy = taxonomy.drop('feature-id', axis=1)
    taxonomy = taxonomy.set_index('NCBI_tax_id')
    taxonomy.index.name = 'Feature ID'

    stratified_table = stratified_table.reset_index()
    stratified_table = stratified_table.drop(columns_to_drop,
                                             axis=1)
    stratified_table = stratified_table.set_index('NCBI_tax_id')
    stratified_table.index.name = 'sample-id'

    return stratified_table.T, taxonomy


plugin.methods.register_function(
    function=select_stratum,
    inputs={'stratified_table': MergedMetaphlanTable},
    parameters={'level': Int % Range(1,None)},
    outputs=[('table', FeatureTable[Frequency]),
             ('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'stratified_table': ('A stratified Metaphlan3 feature table.'),
    },
    parameter_descriptions={
        'level': ('The level (or stratum) of the feature metadata heirarchy '
                  'to select from the input table.')
    },
    output_descriptions={
        'table': ('Filtered table containing only features at specified '
                  'level (or stratum).'),
        'taxonomy': ('Taxonomic feature metadata.')},
    name='Filter Metaphlan3 feature table to single level (or stratum).',
    description=("Filter a Metaphlan3 feature table to the specified "
                 "taxonomic level (or stratum)."),
    citations=[
        citations['bioBakery3']]
)