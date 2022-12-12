# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from qiime2.plugin.testing import TestPluginBase

from q2_sapienns import MetaphlanMergedAbundanceFormat
from q2_sapienns._metaphlan import metaphlan_taxon, frequency


class MetaphlanTaxonTests(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_metaphlan_taxon_bad_level(self):
        _, input_table_df = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            'metaphlan-merged-abundance-1.tsv'
        )

        with self.assertRaisesRegex(ValueError, 'exactly 42 taxonomic levels'):
            metaphlan_taxon(input_table_df, level=42)

    def test_metaphlan_taxon_level_1(self):
        _, input_table_df = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            'metaphlan-merged-abundance-1.tsv')

        obs_table, obs_tax = metaphlan_taxon(input_table_df, level=1)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 2))
        self.assertEqual(list(obs_table.index),
                         ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['k__Archaea',
             'k__Bacteria'])

        self.assertEqual(list(obs_table.T['sample1']),
                         [9.75907, 90.24093])

        self.assertEqual(list(obs_table.T['sample_2']),
                         [0.02352, 99.97648])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (2, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['k__Archaea',
             'k__Bacteria'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['k__Archaea',
             'k__Bacteria'])

    def test_metaphlan_taxon_level_7(self):
        _, input_table_df = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            'metaphlan-merged-abundance-1.tsv')

        obs_table, obs_tax = metaphlan_taxon(input_table_df, level=7)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 7))
        self.assertEqual(list(obs_table.index),
                         ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

        self.assertEqual(list(obs_table.T['sample1']),
                         [9.75907, 45.0, 5.0, 1.24093, 25.0, 10.0, 4.0])

        self.assertEqual(list(obs_table.T['sample_2']),
                         [0.02352, 10.97648, 0.0, 0.0, 0.0, 9.0, 80.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (7, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobrevibacter; s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinobaculum; s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

    def test_metaphlan_taxon_level_7_no_tax_id(self):
        _, input_table_df = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            'metaphlan-merged-abundance-6.tsv')

        obs_table, obs_tax = metaphlan_taxon(input_table_df, level=7)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 7))
        self.assertEqual(list(obs_table.index),
                         ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

        self.assertEqual(list(obs_table.T['sample1']),
                         [9.75907, 45.0, 5.0, 1.24093, 25.0, 10.0, 4.0])

        self.assertEqual(list(obs_table.T['sample_2']),
                         [0.02352, 10.97648, 0.0, 0.0, 0.0, 9.0, 80.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (7, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobrevibacter; s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinobaculum; s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces; s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

    def test_frequency_100000(self):
        _, input_table_df = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            'metaphlan-merged-abundance-1.tsv')

        obs_rf_table, _ = metaphlan_taxon(input_table_df, level=7)
        obs_f_table = frequency(obs_rf_table)

        # Assess resulting tables
        self.assertEqual(obs_f_table.index.name, 'sample-id')
        self.assertEqual(obs_f_table.shape, (2, 7))
        self.assertEqual(list(obs_f_table.index),
                         ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_f_table.columns),
            ['k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

        self.assertEqual(list(obs_f_table.T['sample1']),
                         [9759, 45000, 5000, 1241, 25000, 10000, 4000])

        self.assertEqual(list(obs_f_table.T['sample_2']),
                         [24, 10976, 0, 0, 0, 9000, 80000])

    def test_frequency_1(self):
        _, input_table_df = self.transform_format(
            MetaphlanMergedAbundanceFormat, pd.DataFrame,
            'metaphlan-merged-abundance-1.tsv')

        obs_rf_table, _ = metaphlan_taxon(input_table_df, level=7)
        obs_f_table = frequency(obs_rf_table, 1)

        # Assess resulting tables
        self.assertEqual(obs_f_table.index.name, 'sample-id')
        self.assertEqual(obs_f_table.shape, (2, 7))
        self.assertEqual(list(obs_f_table.index),
                         ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_f_table.columns),
            ['k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|s__Actinobaculum_sp_oral_taxon_183',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_graevenitzii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_naeslundii',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_odontolyticus',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris',  # noqa: E501
             'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_sp_HMSC035G02'])  # noqa: E501

        self.assertEqual(list(obs_f_table.T['sample1']),
                         [0, 0, 0, 0, 0, 0, 0])

        self.assertEqual(list(obs_f_table.T['sample_2']),
                         [0, 0, 0, 0, 0, 0, 1])
