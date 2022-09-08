# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from qiime2.plugin.testing import TestPluginBase

from q2_sapienns import HumannGeneFamilyFormat, HumannPathAbundanceFormat
from q2_sapienns._humann import humann_genefamily, humann_pathway


class HumannTests(TestPluginBase):
    package = 'q2_sapienns.tests'


class HumannGeneFamilyTests(HumannTests):

    def test_humann_genefamilies(self):
        _, input_table_df = self.transform_format(
            HumannGeneFamilyFormat, pd.DataFrame, 'humann-genefamilies-1.tsv'
        )

        obs_table, obs_tax = humann_genefamily(input_table_df)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (1, 8))
        self.assertEqual(list(obs_table.index), ['sample1'])
        self.assertEqual(
            list(obs_table.columns),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris'])  # noqa: E501

        self.assertEqual(list(obs_table.T['sample1']),
                         [150.0, 57.0, 5.0, 4.0,
                          1.0, 31.0, 22.0, 7.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (8, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris'])  # noqa: E501

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UniRef50_unknown; g__Bacteroides; s__Bacteroides_fragilis',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_stercoris'])  # noqa: E501

    def test_humann_genefamilies_unchanged_sample_ids(self):
        _, input_table_df = self.transform_format(
            HumannGeneFamilyFormat, pd.DataFrame, 'humann-genefamilies-1.tsv'
        )

        obs_table, obs_tax = humann_genefamily(
            input_table_df, strip_units_from_sample_ids=False)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (1, 8))
        self.assertEqual(list(obs_table.index), ['sample1_Abundance-RPKs'])
        self.assertEqual(
            list(obs_table.columns),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris'])  # noqa: E501

        self.assertEqual(list(obs_table.T['sample1_Abundance-RPKs']),
                         [150.0, 57.0, 5.0, 4.0,
                          1.0, 31.0, 22.0, 7.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (8, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris'])  # noqa: E501

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UniRef50_unknown; g__Bacteroides; s__Bacteroides_fragilis',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_stercoris'])  # noqa: E501

    def test_humann_genefamilies_multi_sample(self):
        _, input_table_df = self.transform_format(
            HumannGeneFamilyFormat, pd.DataFrame, 'humann-genefamilies-2.tsv'
        )

        obs_table, obs_tax = humann_genefamily(input_table_df)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 10))
        self.assertEqual(list(obs_table.index), ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_unknown|g__Bacteroides.s__Bacteroides_finegoldii',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|unclassified'])

        self.assertEqual(list(obs_table.T['sample1']),
                         [150.0, 25.0, 57.0, 5.0, 4.0,
                          1.0, 31.0, 22.0, 7.0, 0.0])

        self.assertEqual(list(obs_table.T['sample_2']),
                         [0.0, 25.0, 0.0, 0.0, 0.0, 0.0,
                          40.0, 2.0, 2.0, 1.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (10, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_unknown|g__Bacteroides.s__Bacteroides_finegoldii',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|unclassified'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UniRef50_unknown; g__Bacteroides; s__Bacteroides_fragilis',
             'UniRef50_unknown; g__Bacteroides; s__Bacteroides_finegoldii',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; unclassified'])

    def test_humann_genefamilies_multi_sample_unchanged_sample_ids(self):
        _, input_table_df = self.transform_format(
            HumannGeneFamilyFormat, pd.DataFrame, 'humann-genefamilies-2.tsv'
        )

        obs_table, obs_tax = humann_genefamily(
            input_table_df, strip_units_from_sample_ids=False
        )

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 10))
        self.assertEqual(list(obs_table.index),
                         ['sample1_Abundance-RPKs', 'sample_2_Abundance-RPKs'])
        self.assertEqual(
            list(obs_table.columns),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_unknown|g__Bacteroides.s__Bacteroides_finegoldii',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|unclassified'])

        self.assertEqual(list(obs_table.T['sample1_Abundance-RPKs']),
                         [150.0, 25.0, 57.0, 5.0, 4.0,
                          1.0, 31.0, 22.0, 7.0, 0.0])

        self.assertEqual(list(obs_table.T['sample_2_Abundance-RPKs']),
                         [0.0, 25.0, 0.0, 0.0, 0.0, 0.0,
                          40.0, 2.0, 2.0, 1.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (10, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis',
             'UniRef50_unknown|g__Bacteroides.s__Bacteroides_finegoldii',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase|unclassified'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UniRef50_unknown; g__Bacteroides; s__Bacteroides_fragilis',
             'UniRef50_unknown; g__Bacteroides; s__Bacteroides_finegoldii',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_fragilis',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; g__Bacteroides; s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon; unclassified',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_vulgatus',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_thetaiotaomicron',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; g__Bacteroides; s__Bacteroides_stercoris',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase; unclassified'])

    def test_humann_genefamilies_multi_sample_destratify(self):
        _, input_table_df = self.transform_format(
            HumannGeneFamilyFormat, pd.DataFrame, 'humann-genefamilies-2.tsv'
        )

        obs_table, obs_tax = humann_genefamily(input_table_df, destratify=True)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 4))
        self.assertEqual(list(obs_table.index), ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['UNMAPPED',
             'UniRef50_unknown',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase'])

        self.assertEqual(list(obs_table.T['sample1']),
                         [187.0, 150.0, 67.0, 60.0])

        self.assertEqual(list(obs_table.T['sample_2']),
                         [42.0, 50.0, 0.0, 45.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (4, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UNMAPPED',
             'UniRef50_unknown',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase'])
        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UNMAPPED',
             'UniRef50_unknown',
             'UniRef50_A6L0N6: Conserved protein found in conjugate transposon',  # noqa: E501
             'UniRef50_O83668: Fructose-bisphosphate aldolase'])


class HumannPathwayTests(HumannTests):
    def test_humann_pathway(self):
        _, input_table_df = self.transform_format(
            HumannPathAbundanceFormat, pd.DataFrame,
            'humann-pathabundance-1.tsv'
        )

        obs_table, obs_tax = humann_pathway(input_table_df)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (1, 9))
        self.assertEqual(list(obs_table.index), ['sample1'])
        self.assertEqual(
            list(obs_table.columns),
            ['UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
             'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
             'UNINTEGRATED|unclassified',
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation|unclassified',
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|unclassified'])

        self.assertEqual(list(obs_table.T['sample1']),
                         [23.0, 20.0, 12.0, 32.5, 4.5,
                          3.0, 16.7, 8.0, 6.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (9, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
             'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
             'UNINTEGRATED|unclassified',
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation|unclassified',
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|unclassified'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UNINTEGRATED; g__Bacteroides; s__Bacteroides_caccae',
             'UNINTEGRATED; g__Bacteroides; s__Bacteroides_finegoldii',
             'UNINTEGRATED; unclassified',
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation; unclassified',
             'PWY-5484: glycolysis II (from fructose-6P); g__Bacteroides; s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P); g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P); unclassified'])

    def test_humann_pathway_unchanged_sample_ids(self):
        _, input_table_df = self.transform_format(
            HumannPathAbundanceFormat, pd.DataFrame,
            'humann-pathabundance-1.tsv'
        )

        obs_table, obs_tax = humann_pathway(
            input_table_df, strip_units_from_sample_ids=False)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (1, 9))
        self.assertEqual(list(obs_table.index), ['sample1_Abundance'])
        self.assertEqual(
            list(obs_table.columns),
            ['UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
             'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
             'UNINTEGRATED|unclassified',
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation|unclassified',
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|unclassified'])

        self.assertEqual(list(obs_table.T['sample1_Abundance']),
                         [23.0, 20.0, 12.0, 32.5, 4.5,
                          3.0, 16.7, 8.0, 6.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (9, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
             'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
             'UNINTEGRATED|unclassified',
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation|unclassified',
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|unclassified'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UNINTEGRATED; g__Bacteroides; s__Bacteroides_caccae',
             'UNINTEGRATED; g__Bacteroides; s__Bacteroides_finegoldii',
             'UNINTEGRATED; unclassified',
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation; unclassified',
             'PWY-5484: glycolysis II (from fructose-6P); g__Bacteroides; s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P); g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P); unclassified'])

    def test_humann_pathway_multisample(self):
        _, input_table_df = self.transform_format(
            HumannPathAbundanceFormat, pd.DataFrame,
            'humann-pathabundance-2.tsv'
        )

        obs_table, obs_tax = humann_pathway(input_table_df)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 10))
        self.assertEqual(list(obs_table.index), ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
             'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
             'UNINTEGRATED|unclassified',
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'PWY0-1301: melibiose degradation|unclassified',
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|unclassified'])

        self.assertEqual(list(obs_table.T['sample1']),
                         [23.0, 20.0, 12.0, 32.5, 4.5,
                          0.0, 3.0, 16.7, 8.0, 6.0])
        self.assertEqual(list(obs_table.T['sample_2']),
                         [0.0, 23.0, 8.0, 10.0, 0.0, 1.0,
                          11.0, 5.0, 2.0, 4.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (10, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae',
             'UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii',
             'UNINTEGRATED|unclassified',
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_fragilis',  # noqa: E501
             'PWY0-1301: melibiose degradation|unclassified',
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P)|unclassified'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UNINTEGRATED; g__Bacteroides; s__Bacteroides_caccae',
             'UNINTEGRATED; g__Bacteroides; s__Bacteroides_finegoldii',
             'UNINTEGRATED; unclassified',
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_caccae',  # noqa: E501
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'PWY0-1301: melibiose degradation; g__Bacteroides; s__Bacteroides_fragilis',  # noqa: E501
             'PWY0-1301: melibiose degradation; unclassified',
             'PWY-5484: glycolysis II (from fructose-6P); g__Bacteroides; s__Bacteroides_caccae',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P); g__Bacteroides; s__Bacteroides_finegoldii',  # noqa: E501
             'PWY-5484: glycolysis II (from fructose-6P); unclassified'])

    def test_humann_pathway_multisample_destratify(self):
        _, input_table_df = self.transform_format(
            HumannPathAbundanceFormat, pd.DataFrame,
            'humann-pathabundance-2.tsv'
        )

        obs_table, obs_tax = humann_pathway(input_table_df, destratify=True)

        # Assess resulting tables
        self.assertEqual(obs_table.index.name, 'sample-id')
        self.assertEqual(obs_table.shape, (2, 4))
        self.assertEqual(list(obs_table.index), ['sample1', 'sample_2'])
        self.assertEqual(
            list(obs_table.columns),
            ['UNMAPPED',
             'UNINTEGRATED',
             'PWY0-1301: melibiose degradation',
             'PWY-5484: glycolysis II (from fructose-6P)'])

        self.assertEqual(list(obs_table.T['sample1']),
                         [140.0, 87.0, 57.5, 54.7])
        self.assertEqual(list(obs_table.T['sample_2']),
                         [99.0, 42.0, 42.0, 99.0])

        # Assess resulting feature metadata
        self.assertEqual(obs_tax.index.name, 'Feature ID')
        self.assertEqual(obs_tax.shape, (4, 1))
        self.assertEqual(list(obs_tax.columns), ['Taxon'])
        self.assertEqual(
            list(obs_tax.index),
            ['UNMAPPED',
             'UNINTEGRATED',
             'PWY0-1301: melibiose degradation',
             'PWY-5484: glycolysis II (from fructose-6P)'])

        self.assertEqual(
            list(obs_tax['Taxon']),
            ['UNMAPPED',
             'UNINTEGRATED',
             'PWY0-1301: melibiose degradation',
             'PWY-5484: glycolysis II (from fructose-6P)'])
