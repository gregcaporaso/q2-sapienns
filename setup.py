# ----------------------------------------------------------------------------
# Copyright (c) 2022, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name='q2-sapienns',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Greg Caporaso",
    author_email="greg.caporaso@nau.edu",
    description=(
        'A plugin for interacting with data generated by the biobakery '
        'tools.'
    ),
    url="https://github.com/q2-sapienns",
    entry_points={
        'qiime2.plugins':
        ['q2-sapienns=q2_sapienns.plugin_setup:plugin']
    },
    package_data={
        'q2_sapienns': ['citations.bib'],
        'q2_sapienns.tests': ['data/*', 'data/*/*'],
        'q2_sapienns.types.tests': ['data/*']
    },
    zip_safe=False,
)