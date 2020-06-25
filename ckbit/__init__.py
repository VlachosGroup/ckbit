# -*- coding: utf-8 -*-

####
#
# setuptools likes to see a name for the package,
# and it's best-practices to have the __version__
# present, too:
#
name = 'ckbit'
__version__ = '1.0.0'

# Pull rxn_ord, pfr, cstr, app_ea in as ckbit.<x>:
from . import rxn_ord
from . import pfr
from . import cstr
from . import app_ea
