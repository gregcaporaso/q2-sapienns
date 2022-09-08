# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from q2_sapienns.plugin_setup import plugin as sapienns_plugin


class PluginSetupTests(unittest.TestCase):

    def test_plugin_setup(self):
        self.assertEqual(sapienns_plugin.name, 'sapienns')
