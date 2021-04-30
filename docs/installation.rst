.. _installation:

Installation
************

Installation Options Overview
-----------------------------
The user has two different installation options available:

1. Through the use of a pip command

2. Through the use of a Docker container image

We recommend most users proceed with option 1 so they may interface with the 
Python code as they would standardly. However, for users who are experienced
with containers, the container provides a more simplified installation route.
CKBIT requires more steps to install than a simple pip command due to the 
complexity of the dependency providing the Bayesian inference capabilities.
We guide the user through these steps below.

--------------------------------------------------------------------------------

Installation Option 1
---------------------

**Steps Overview**

1. Installing Python
2. Installing PyStan
3. Installing CKBIT

**Installing Python**

Anaconda is the recommended method to install Python for scientific 
applications. It is supported on Linux, Windows and Mac OS X. This software 
runs on Python 3.7 and more recent versions. Anaconda can be downloaded here:
https://www.anaconda.com/products/individual#download-section

**Installing PyStan**

PyStan is the Python library that provides the Bayesian inference functionality
for CKBIT, so it is crucial to install it properly. This becomes complicated
since there are multiple versions of PyStan with different installation pathways
depending on whether the computer is running Linux, Windows, or a Mac OS. To run 
CKBIT, **a user must first install PyStan with the CVODES library. Without the
CVODES version of PyStan, CKBIT will not install.** For Windows users, please 
install Anaconda and see below for instructions. For all other users, please 
reference the documentation provided by PyStan on how to install their 
software here: 
https://pystan2.readthedocs.io/en/latest/installation_beginner.html

For details on how to get the CVODES version, reference this link:
https://pystan2.readthedocs.io/en/latest/installation_cvodes.html

For Windows users with Anaconda, follow these steps:

1. Open the Anaconda prompt application.

2. Ensure you have the necessary dependencies for PyStan: m2w64-toolchain,
numpy, cython, and git. To do this, enter 
::

	conda list
	
into the Anaconda command prompt. The output will be a list of all dependences 
in the current environment. Ensure that all the listed dependencies are present. 
If any are missing, see below. It is also important to confirm that the cython 
version is >=0.22 and the numpy version is >=1.7. If the versions
of these dependencies are out of date, update them by entering the following 
into the Anaconda command prompt (replacing dependency_name with the name of 
the dependency you wish to update):
::

	conda update dependency_name
	
If any of the dependencies are missing from the list, they simply need to be 
installed. The following commands entered into the Anaconda command prompt can 
be used for each dependency:

m2w64-toolchain:
::

	conda install libpython m2w64-toolchain -c msys2

numpy:
::

	conda install -c anaconda numpy 

cython:
::

	conda install -c anaconda cython 
	 
git:
::

	conda install -c anaconda git 

3. The final step is to install the CVODES version of PyStan. Enter the 
following commands into the Anaconda command prompt sequentially hitting 
y for yes when appropriate:
::

	git config --system core.longpaths true
	git clone --recursive https://github.com/stan-dev/pystan2 
	cd pystan2
	git checkout cvodes
	git submodule update --recursive
	python -c "import os, shutil; [shutil.rmtree(r'\\?\{}'.format(os.path.abspath(dirname)), ignore_errors=True) for dirname in [dirname for dirname, *_ in os.walk('pystan/stan') if any(dirname.endswith(ends) for ends in ['doc', 'test'])]]"
	python -m pip install .

**Installing CKBIT**

Now that the CVODES version of PyStan is installed, the user can now install
CKBIT. Run the command in the anaconda prompt of:
::

	pip install ckbit

--------------------------------------------------------------------------------

Installation Option 2
---------------------
The container is publicly available on Docker Hub. 

A Docker container image can be accessed at this link:



Once you have successfully downloaded the container:

1. Load the image:
::
	
	docker image load --input ckbit-container-latest.tar

2. Start a Python shell inside a container based off of that image:
::

    docker run --rm --interactive --tty ckbit:0.0.1 python3

3. Interface with the Python prompt to run CKBIT
