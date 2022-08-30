import os
import os.path
import shutil
import unittest

from q2_sapienns import (
    HumannGeneFamilyFormat,
    HumannPathAbundanceFormat,
    MetaphlanMergedAbundanceFormat,
)
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError


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


class TestMetaphlanMergedAbundanceFormat(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_metaphlan_merged_abundance_format_valid(self):
        filenames = ['metaphlan-merged-abundance-1.tsv']
        filepaths = [self.get_data_path(filename)
                     for filename in filenames]

        for filepath in filepaths:
            format = MetaphlanMergedAbundanceFormat(filepath, mode='r')
            format.validate()
