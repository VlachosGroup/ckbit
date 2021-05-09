#!/usr/bin/env python
#
# Setuptools install script for CKBIT
#

import sys
import os
import subprocess
import setuptools

#
# Don't do the extensive checks for PyStan unless we need to:
#
if not( len(sys.argv) >= 2 and ('--help' in sys.argv[1:] or sys.argv[1] in ('--help-commands', 'egg_info', '--version', 'clean')) ):
    #
    # Test for presence of the PyStan library.  Note that we DO NOT exit if
    # the tests do not succeed, we merely remind the user that PyStan with
    # CVODES is a requirement s/he will need to satisfy.
    #
    try:
        import pystan
        try:
            #
            # So PyStan is present...does it have CVODES support?  We'll
            # run that test as a forked process so we can capture stderr/stdout
            # and probably not display it -- the PyStan compilation can be a
            # tad verbose!
            #
            sp = subprocess.Popen(
                        args=(sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pystan-test.py')),
                        stdin=None,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
            (sp_out, sp_err) = sp.communicate()
            if sp.returncode != 0:
                #
                # So pystan is present but running a CVODES model failed.
                #
                sys.stderr.write("""
**
** WARNING:  this software requires a PyStan package built with CVODES
**           functionality.  Copies of PyStan on most Conda channels do not
**           include this functionality, so pystan is not a dependency that
**           can be satisfied by pip or conda installs.
**
**           If the accessible PyStan package does not include CVODES
**           functionality you will see runtime exceptions when Stan code
**           is compiled.  Please see the pystan-install.sh script for a
**           working method for building a copy of PyStan
**

""")
        except Exception as E:
            sys.stderr.write('ERROR:  failed to execute PyStan CVODES test: {:s}\n'.format(str(E)))
            sys.exit(1)
    except:
        sys.stderr.write("""
**
** WARNING:  this software requires a PyStan package built with CVODES
**           functionality.  Please see the pystan-install.sh script
**           for a working method for building a copy of PyStan.  Copies of
**           PyStan on most Conda channels do not include this functionality,
**           so PyStan is not a dependency that can be satisfied by pip or
**           conda installs.
**

""")

setuptools_info = {
    'name': 'ckbit',
    'version': '1.0.0',
    'author': 'Vlachos Research Group',
    'author_email': 'vlachos@udel.edu',
    'description': 'Kinetic Bayesian Inference',
    'zip_safe': True,
    'url': 'https://github.com/VlachosGroup/ckbit',
    'packages': setuptools.find_packages(),
    'python_requires': '>=3.7',
    'install_requires': [
        'numpy>=1.16',
        'datetime',
        'tabulate>=0.8',
        'VUnits>=0.0.3',
        'arviz>=0.4',
        'pandas>=0.25',
        'xlrd',
        'openpyxl>=3.0.0',
        'seaborn>=0.9.0',
        'matplotlib>=3.1',
        ],
    'classifiers': [
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        ],
    }
if os.path.isfile('README.rst'):
    with open('README.rst', 'r') as fh:
        setuptools_info['long_description'] = fh.read()
        if sys.version_info[0] >= 3:
            #
            # Augment for Python 3 setuptools:
            #
            setuptools_info['long_description_content_type'] = 'text/x-rst'

setuptools.setup(**setuptools_info)
