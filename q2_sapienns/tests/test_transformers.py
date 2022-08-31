import pandas as pd
from pandas.testing import assert_frame_equal

from qiime2.plugin.testing import TestPluginBase

from q2_sapienns import (
    HumannGeneFamilyFormat, HumannPathAbundanceFormat,
    MetaphlanMergedAbundanceFormat)


class TestHumannFormatTransformers(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_gene_family_format_to_dataframe_with_header(self):
        index = ['UNMAPPED',
                 'UniRef50_unknown',
                 'UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
                 'UniRef50_unknown|g__Bacteroides.s__Bacteroides_finegoldii',
                 'UniRef50_A6L0N6: Conserved protein found in conjugate transposon',  # noqa: E501
                 'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
                 'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
                 'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
                 'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
                 'UniRef50_O83668: Fructose-bisphosphate aldolase',
                 'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
                 'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
                 'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
                 'UniRef50_O83668: Fructose-bisphosphate aldolase|unclassified']  # noqa: E501
        index = pd.Index(index, name='feature-id', dtype=object)
        columns = ['sample1_Abundance-RPKs', 'sample_2_Abundance-RPKs']
        exp = pd.DataFrame([[187.0, 42.0],
                            [150.0, 50.0],
                            [150.0, 0.0],
                            [25.0, 25.0],
                            [67.0, 0.0],
                            [57.0, 0.0],
                            [5.0, 0.0],
                            [4.0, 0.0],
                            [1.0, 0.0],
                            [60.0, 45.0],
                            [31.0, 40.0],
                            [22.0, 2.0],
                            [7.0, 2.0],
                            [0.0, 1.0]],
                           index=index,
                           columns=columns,
                           dtype=float)

        _, obs = self.transform_format(
            HumannGeneFamilyFormat, pd.DataFrame,
            filename='humann-genefamilies-2.tsv')

        assert_frame_equal(obs, exp)

    def test_path_abundance_format_to_dataframe_with_header(self):
        index = ['UNMAPPED',
                 'UNINTEGRATED',
                 'UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
                 'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
                 'UNINTEGRATED|unclassified',
                 'PWY0-1301: melibiose degradation',
                 'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
                 'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
                 'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
                 'PWY0-1301: melibiose degradation|unclassified',
                 'PWY-5484: glycolysis II (from fructose-6P)',
                 'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
                 'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
                 'PWY-5484: glycolysis II (from fructose-6P)|unclassified']
        index = pd.Index(index, name='feature-id', dtype=object)
        columns = ['sample1_Abundance', 'sample_2_Abundance']
        exp = pd.DataFrame([[140.0, 99.0],
                            [87.0, 42.0],
                            [23.0, 0.0],
                            [20.0, 23.0],
                            [12.0, 8.0],
                            [57.5, 42.0],
                            [32.5, 10.0],
                            [4.5, 0.0],
                            [0.0, 1.0],
                            [3.0, 11.0],
                            [54.7, 99.0],
                            [16.7, 5.0],
                            [8.0, 2.0],
                            [6.0, 4.0]],
                           index=index,
                           columns=columns,
                           dtype=float)

        _, obs = self.transform_format(
            HumannPathAbundanceFormat, pd.DataFrame,
            filename='humann-pathabundance-2.tsv')

        assert_frame_equal(obs, exp)


class TestMetaphlanFormatTransformers(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_metaphlan_format_to_dataframe_with_header(self):
        index = ['k__Archaea',
                 'k__Archaea|p__Euryarchaeota',
                 'k__Archaea|p__Euryarchaeota|c__Methanobacteria',
                 'k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales',  # noqa: E501
                 'k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae',  # noqa: E501
                 'k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter',  # noqa: E501
                 'k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
                 'k__Bacteria',
                 'k__Bacteria|p__Actinobacteria',
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria',
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
                 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02']  # noqa: E501
        index = pd.Index(index, name='feature-id', dtype=object)
        columns = ['NCBI_tax_id', 'sample1', 'sample_2']
        exp = pd.DataFrame([['2157', 9.75907, 0.02352],
                            ['2157|28890', 9.75907, 0.02352],
                            ['2157|28890|183925', 9.75907, 0.02352],
                            ['2157|28890|183925|2158', 9.75907, 0.02352],
                            ['2157|28890|183925|2158|2159', 9.75907, 0.02352],
                            ['2157|28890|183925|2158|2159|2172', 9.75907, 0.02352],  # noqa: E501
                            ['2157|28890|183925|2158|2159|2172|2173', 9.75907, 0.02352],  # noqa: E501
                            ['2', 90.24093, 99.97648],
                            ['2|201174', 90.24093, 99.97648],
                            ['2|201174|1760', 90.24093, 99.97648],
                            ['2|201174|1760|2037', 90.24093, 99.97648],
                            ['2|201174|1760|2037|2049', 90.24093, 99.97648],
                            ['2|201174|1760|2037|2049|76833', 45.0, 10.97648],
                            ['2|201174|1760|2037|2049|76833|712888', 45.0, 10.97648],  # noqa: E501
                            ['2|201174|1760|2037|2049|1654', 45.24093, 89],
                            ['2|201174|1760|2037|2049|1654|55565', 5.0, 0.0],
                            ['2|201174|1760|2037|2049|1654|1655', 1.24093, 0.0],  # noqa: E501
                            ['2|201174|1760|2037|2049|1654|1660', 25.0, 0.0],
                            ['2|201174|1760|2037|2049|1654|544580', 10.0, 9.0],
                            ['2|201174|1760|2037|2049|1654|1739406', 4.0, 80.0]],  # noqa: E501
                           index=index,
                           columns=columns)

        _, obs = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            filename='metaphlan-merged-abundance-1.tsv')

        assert_frame_equal(obs, exp)
