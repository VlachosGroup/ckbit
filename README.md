# KBI

Kinetic Bayesian Inference

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Installation

The `setup.py` file is a setuptools installation script.  It can be used to install into a Python virtual environment, a specific directory, or right into your native Python install.

The CKBIT module has several external requirements, most of which can be satisfied by Conda or PIP:

- python (3.7 or newer)
- numpy
- cython (0.22 or newer)
- datetime
- tabulate
- VUnits
- arviz
- pandas
- xlrd
- matplotlib
- PyStan (**with CVODES functionality**)

The final dependency is the most difficult:  while many conda channels contain PyStan packages, none include the additional CVODES functionality required by CKBIT.  A shell script (`pystan-install.sh`) is included that was used to build and install a CVODES-enabled version of PyStan.

### Virtual Environment

A fresh Python virtual environment can be created using Intel Python or Anaconda:

```
$ conda create --prefix="${HOME}/pystan-env" python'>=3.7' pip numpy cython'>=0.22'
   :
$ source activate "${HOME}/pystan-env"
```

Enter the CKBIT source directory and execute the PyStan+CVODES build:

```
(pystan-env)$ cd CKBIT-1.0.0
(pystan-env)$ ./pystan-install.sh --verbose
OK - cloned git repository
OK - checked-out cvodes branch
OK - updated submodules in repository
OK - removed unneeded components from repository
OK - PyStan+CVODES module built
OK - PyStan+CVODES module installed
```

Note that the GNU C++ compiler will require a **lot** of memory to process the Stan source code; in a Docker container (`continuumio/anaconda3:latest`) upwards of 9 GiB of RAM was required for the build to succeed.

At this point if the installation was successful the module can be tested:

```
(pystan-env)$ python pystan-test.py 
    :
(pystan-env)$ echo $?
0
```

You are likely to see a few compiler warnings as the Stan code is generated and compiled.  Ideally that is all you should see and the script return code should be zero.  If so, then you're ready to install CKBIT into the same Python virtual environment:

```
(pystan-env)$ python setup.py install
   :
Using {..}pystan-env/lib/python3.8/site-packages
Finished processing dependencies for ckbit==1.0.0
```

This will take a while at first since it, too, is going to run the `python-test.py` script to see if you have an appropriate PyStan+CVODES present in the environment into which you are installing.  If not found, the script will display a warning about getting one installed into the environment — but it does not prevent the installation of CKBIT or its other dependencies into the environment.

You should now have a usable CKBIT virtual environment:

```
(pystan-env)$ python -c 'import ckbit; print(ckbit);'
<module 'ckbit' from '{..}/CKBIT-1.0.0/ckbit/__init__.py'>
```
