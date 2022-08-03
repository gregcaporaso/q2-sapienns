from qiime2.plugin import Plugin

import q2_sapienns

plugin = Plugin(
    name='sapienns',
    version=q2_sapienns.__version__,
    website='https://qiime2.org',
    user_support_text='https://forum.qiime2.org',
    package='q2_sapienns'
)


