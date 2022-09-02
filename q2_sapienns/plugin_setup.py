from qiime2.plugin import (Plugin, SemanticType, TextFileFormat, model,
                           ValidationError, Citations, Int, Range, Bool)
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.feature_data import FeatureData, Taxonomy

import q2_sapienns
from ._humann import humann_pathway, humann_genefamily
from ._metaphlan import metaphlan_taxon

import pandas as pd

plugin = Plugin(
    name='sapienns',
    version=q2_sapienns.__version__,
    website='https://qiime2.org',
    user_support_text='https://forum.qiime2.org',
    package='q2_sapienns'
)

MetaphlanMergedAbundanceTable = SemanticType('MetaphlanMergedAbundanceTable')
HumannPathAbundanceTable = SemanticType('HumannPathAbundanceTable')
HumannGeneFamilyTable = SemanticType('HumannGeneFamilyTable')

plugin.register_semantic_types(MetaphlanMergedAbundanceTable)
plugin.register_semantic_types(HumannPathAbundanceTable)
plugin.register_semantic_types(HumannGeneFamilyTable)


class MetaphlanMergedAbundanceFormat(TextFileFormat):

    def _equal_number_of_columns(self, n_lines):
        with self.open() as fh:
            header_line = fh.readline()
            while header_line.startswith('#'):
                header_line = fh.readline()
            n_header_fields = len(header_line.split('\t'))
            if n_header_fields < 3:
                raise ValidationError(
                    'No sample columns appear to be present.')
            for idx, line in enumerate(fh, 2):
                if n_lines is not None and idx > n_lines + 1:
                    break
                fields = line.strip().split('\t')
                n_fields = len(fields)
                if n_fields != n_header_fields:
                    raise ValidationError(
                        'Number of columns on line %d is inconsistent with '
                        'the header line.' % line)
                for value in fields[2:]:
                    try:
                        value = float(value)
                    except ValueError:
                        raise ValidationError(
                            'Values in table must be float-able. Found: %s' %
                            value
                        )
                    if value > 100.0 or value < 0.0:
                        raise ValidationError(
                            'Values must be in range [0, 100]. Found: %f' %
                            value
                        )

    def _validate_(self, level):
        level_to_n_lines = {'min': 5, 'max': None}
        self._equal_number_of_columns(level_to_n_lines[level])


class HumannTableFormat(TextFileFormat):

    def _equal_number_of_columns(self, n_lines):
        with self.open() as fh:
            header_line = fh.readline().strip()
            header_fields = header_line.split('\t')
            n_header_fields = len(header_fields)
            if n_header_fields < 2:
                raise ValidationError(
                    'No sample columns appear to be present.')
            for sample_id in header_fields[1:]:
                if not sample_id.endswith(self._unit_label):
                    raise ValidationError(
                        'Expected sample ids (e.g., %s) to end with unit '
                        'descriptor %s' % (sample_id, self._unit_label)
                    )
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


class HumannPathAbundanceFormat(HumannTableFormat):
    _unit_label = 'Abundance'


class HumannGeneFamilyFormat(HumannTableFormat):
    _unit_label = 'RPKs'


MetaphlanMergedAbundanceDirectoryFormat = model.SingleFileDirectoryFormat(
    'MetaphlanMergedAbundanceDirectoryFormat', 'table.tsv',
    MetaphlanMergedAbundanceFormat)

HumannPathAbundanceDirectoryFormat = model.SingleFileDirectoryFormat(
    'HumannPathAbundanceDirectoryFormat', 'table.tsv',
    HumannPathAbundanceFormat)

HumannGeneFamilyDirectoryFormat = model.SingleFileDirectoryFormat(
    'HumannGeneFamilyDirectoryFormat', 'table.tsv', HumannGeneFamilyFormat)

plugin.register_formats(MetaphlanMergedAbundanceFormat,
                        MetaphlanMergedAbundanceDirectoryFormat)

plugin.register_semantic_type_to_format(
    MetaphlanMergedAbundanceTable, MetaphlanMergedAbundanceDirectoryFormat)


plugin.register_formats(HumannPathAbundanceFormat,
                        HumannPathAbundanceDirectoryFormat)

plugin.register_semantic_type_to_format(HumannPathAbundanceTable,
                                        HumannPathAbundanceDirectoryFormat)


plugin.register_formats(HumannGeneFamilyFormat,
                        HumannGeneFamilyDirectoryFormat)

plugin.register_semantic_type_to_format(HumannGeneFamilyTable,
                                        HumannGeneFamilyDirectoryFormat)


def _humann_to_df(ff):
    result = pd.read_csv(str(ff), sep='\t', header=0, index_col=0)
    result.index.name = 'feature-id'
    return result


@plugin.register_transformer
def _1(ff: MetaphlanMergedAbundanceFormat) -> pd.DataFrame:
    result = pd.read_csv(str(ff), sep='\t', header=0, index_col=0,
                         comment='#')
    result.index.name = 'feature-id'
    return result


@plugin.register_transformer
def _2(ff: HumannPathAbundanceFormat) -> pd.DataFrame:
    return _humann_to_df(ff)


@plugin.register_transformer
def _3(ff: HumannGeneFamilyFormat) -> pd.DataFrame:
    return _humann_to_df(ff)


citations = Citations.load('citations.bib', package='q2_sapienns')


plugin.methods.register_function(
    function=metaphlan_taxon,
    inputs={'stratified_table': MetaphlanMergedAbundanceTable},
    parameters={'level': Int % Range(1, None)},
    outputs=[('table', FeatureTable[RelativeFrequency]),
             ('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'stratified_table': ('A stratified MetaPhlAn3 feature table.'),
    },
    parameter_descriptions={
        'level': ('The level (or stratum) of the feature metadata heirarchy '
                  'to select from the input table.')
    },
    output_descriptions={
        'table': ('Filtered table containing only features at specified '
                  'level (or stratum).'),
        'taxonomy': ('Taxonomic feature metadata.')},
    name='Filter MetaPhlAn3 feature table to single level (or stratum).',
    description=('Filter a MetaPhlAn3 feature table to the specified '
                 'taxonomic level (or stratum).'),
    citations=[
        citations['bioBakery3']]
)


plugin.methods.register_function(
    function=humann_pathway,
    inputs={'pathway_table': HumannPathAbundanceTable},
    parameters={'strip_units_from_sample_ids': Bool,
                'destratify': Bool},
    outputs=[('table', FeatureTable[Frequency]),
             ('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'pathway_table': ('A stratified HUMAnN3 pathway table.'),
    },
    parameter_descriptions={
        'strip_units_from_sample_ids': 'Remove units from input sample ids.',
        'destratify': ('Only include un-stratified pathways (i.e., those not '
                       'including taxa) in the output table. By default, only '
                       'stratified pathways will be included in the output '
                       'table.')
    },
    output_descriptions={
        'table': ('Output feature table.'),
        'taxonomy': ('Output feature metadata.')},
    name='Prepare HUMAnN3 pathway data.',
    description=('Prepare HUMAnN3 pathway table and pathway metadata for '
                 'QIIME 2. The table is prepared so that only '
                 'stratified (i.e., including taxa) or unstratified (i.e., '
                 'not including taxa) data are included in the output.'),
    citations=[
        citations['bioBakery3']]
)

plugin.methods.register_function(
    function=humann_genefamily,
    inputs={'genefamily_table': HumannGeneFamilyTable},
    parameters={'strip_units_from_sample_ids': Bool,
                'destratify': Bool},
    outputs=[('table', FeatureTable[Frequency]),
             ('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'genefamily_table': ('A stratified HUMAnN3 gene family table.'),
    },
    parameter_descriptions={
        'strip_units_from_sample_ids': 'Remove units from input sample ids.',
        'destratify': ('Only include un-stratified gene families (i.e., those '
                       'not including taxa) in the output table. By default, '
                       'only stratified pathways will be included in the '
                       'output table.')
    },
    output_descriptions={
        'table': ('Output feature table.'),
        'taxonomy': ('Output feature metadata.')},
    name='Prepare HUMAnN3 gene family data.',
    description=('Prepare HUMAnN3 gene family table and gene family metadata '
                 'for QIIME 2. The table is prepared so that only '
                 'stratified (i.e., including taxa) or unstratified (i.e., '
                 'not including taxa) data are included in the output.'),
    citations=[
        citations['bioBakery3']]
)
