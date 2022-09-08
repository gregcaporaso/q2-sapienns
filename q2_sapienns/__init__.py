# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# flake8: noqa

from . import _version
__version__ = _version.get_versions()['version']

from .plugin_setup import (
    HumannGeneFamilyDirectoryFormat, HumannGeneFamilyFormat,
    HumannPathAbundanceDirectoryFormat, HumannPathAbundanceFormat,
    MetaphlanMergedAbundanceDirectoryFormat, MetaphlanMergedAbundanceFormat,
    MetaphlanMergedAbundanceTable, HumannPathAbundanceTable,
    HumannGeneFamilyTable
)
