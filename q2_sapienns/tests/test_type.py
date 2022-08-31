from q2_sapienns import (
    MetaphlanMergedAbundanceTable, HumannPathAbundanceTable,
    HumannGeneFamilyTable, HumannGeneFamilyDirectoryFormat,
    HumannPathAbundanceDirectoryFormat,
    MetaphlanMergedAbundanceDirectoryFormat
)

from qiime2.plugin.testing import TestPluginBase


class TestTypes(TestPluginBase):
    package = 'q2_sapienns.tests'

    def test_metaphlan_semantic_types_registration(self):
        self.assertRegisteredSemanticType(MetaphlanMergedAbundanceTable)
        self.assertSemanticTypeRegisteredToFormat(
            MetaphlanMergedAbundanceTable,
            MetaphlanMergedAbundanceDirectoryFormat)

    def test_humann_semantic_types_registration(self):
        self.assertRegisteredSemanticType(HumannGeneFamilyTable)
        self.assertSemanticTypeRegisteredToFormat(
            HumannGeneFamilyTable,
            HumannGeneFamilyDirectoryFormat)

        self.assertRegisteredSemanticType(HumannPathAbundanceTable)
        self.assertSemanticTypeRegisteredToFormat(
            HumannPathAbundanceTable,
            HumannPathAbundanceDirectoryFormat)
