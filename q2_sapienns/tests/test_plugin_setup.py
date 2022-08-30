import unittest

from q2_sapienns.plugin_setup import plugin as sapienns_plugin


class PluginSetupTests(unittest.TestCase):

    def test_plugin_setup(self):
        self.assertEqual(sapienns_plugin.name, 'sapienns')