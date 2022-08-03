from qiime2.plugin import (Plugin, SemanticType, TextFileFormat, model,
                           ValidationError, Citations, Int, Range)
from q2_types.feature_table import FeatureTable, Frequency

import q2_sapienns

import pandas as pd

plugin = Plugin(
    name='sapienns',
    version=q2_sapienns.__version__,
    website='https://qiime2.org',
    user_support_text='https://forum.qiime2.org',
    package='q2_sapienns'
)

Bb3StratifiedTable = SemanticType('Bb3StratifiedTable')

plugin.register_semantic_types(Bb3StratifiedTable)

class Bb3StatifiedTableFormat(TextFileFormat):

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

Bb3StatifiedTableDirectoryFormat = model.SingleFileDirectoryFormat(
    'Bb3StatifiedTableDirectoryFormat', 'table.tsv', Bb3StatifiedTableFormat)

plugin.register_formats(Bb3StatifiedTableFormat,
                        Bb3StatifiedTableDirectoryFormat)

plugin.register_semantic_type_to_format(Bb3StratifiedTable,
                                        Bb3StatifiedTableDirectoryFormat)


@plugin.register_transformer
def _1(ff: Bb3StatifiedTableFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', header=0, index_col=0)
    return result

citations = Citations.load('citations.bib', package='q2_sapienns')

def table_at_level(stratified_table: pd.DataFrame, level: int) -> pd.DataFrame:
    # TODO: this doesn't actually do anything interest yet
    return stratified_table.T


plugin.methods.register_function(
    function=table_at_level,
    inputs={'stratified_table': Bb3StratifiedTable},
    parameters={'level': Int % Range(1,None)},
    outputs=[('table', FeatureTable[Frequency])],
    input_descriptions={
        'stratified_table': ('A stratified bioBakery 3 feature table.'),
    },
    parameter_descriptions={
        'level': ('The level of the feature metadata heirarchy to select from'
                  ' the input table.')
    },
    output_descriptions={'table':
     ('Filtered table containing only features at specified level.')},
    name='Filter bioBakery3 feature table to single level.',
    description=("Filter a bioBakery 3 feature table to the specified level."),
    citations=[
        citations['bioBakery3']]
)