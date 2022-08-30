
from . import _version
__version__ = _version.get_versions()['version']

from .plugin_setup import (HumannGeneFamilyDirectoryFormat,
                           HumannGeneFamilyFormat,
)
