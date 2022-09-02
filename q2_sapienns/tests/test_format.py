from q2_sapienns import (
    HumannGeneFamilyFormat,
    HumannPathAbundanceFormat,
    MetaphlanMergedAbundanceFormat,
)
from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase


class TestHumannGeneFamilyFormat(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_genefamily_format_valid(self):
        filenames = ['humann-genefamilies-1.tsv',
                     'humann-genefamilies-2.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            format = HumannGeneFamilyFormat(filepath, mode='r')
            format.validate()

    def test_genefamily_format_invalid_ids(self):
        filenames = ['humann-genefamilies-3.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, r's1\) .* RPKs'):
                format = HumannGeneFamilyFormat(filepath, mode='r')
                format.validate()

    def test_genefamily_format_no_samples(self):
        filenames = ['humann-genefamilies-4.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, 'No sample columns'):
                format = HumannGeneFamilyFormat(filepath, mode='r')
                format.validate()


class TestHumannPathAbundanceFormat(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_pathabundance_format_valid(self):
        filenames = ['humann-pathabundance-1.tsv',
                     'humann-pathabundance-2.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            format = HumannPathAbundanceFormat(filepath, mode='r')
            format.validate()

    def test_pathabundance_format_invalid_ids(self):
        filenames = ['humann-pathabundance-3.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, r's1\) .* Abundance'):
                format = HumannPathAbundanceFormat(filepath, mode='r')
                format.validate()

    def test_pathabundance_format_no_samples(self):
        filenames = ['humann-pathabundance-4.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, 'No sample columns'):
                format = HumannPathAbundanceFormat(filepath, mode='r')
                format.validate()


class TestMetaphlanMergedAbundanceFormat(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_metaphlan_merged_abundance_format_valid(self):
        filenames = ['metaphlan-merged-abundance-1.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            format = MetaphlanMergedAbundanceFormat(filepath, mode='r')
            format.validate()

    def test_metaphlan_merged_abundance_format_no_samples(self):
        filenames = ['metaphlan-merged-abundance-2.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, 'No sample columns'):
                format = MetaphlanMergedAbundanceFormat(filepath, mode='r')
                format.validate()

    def test_metaphlan_merged_abundance_format_value_out_of_range(self):
        filenames = ['metaphlan-merged-abundance-3.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, r'range .* 100\.001'):
                format = MetaphlanMergedAbundanceFormat(filepath, mode='r')
                format.validate()

    def test_metaphlan_merged_abundance_format_invalid_value_type(self):
        filenames = ['metaphlan-merged-abundance-4.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            with self.assertRaisesRegex(ValidationError, 'float.*gigawatt'):
                format = MetaphlanMergedAbundanceFormat(filepath, mode='r')
                format.validate()
