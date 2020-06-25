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
import rxn_ord
import pfr
import cstr
import app_ea
