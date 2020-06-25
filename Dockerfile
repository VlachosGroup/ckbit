#
# This is the definition of the CKBIT Docker container.  Based off
# the official Python 3 Anaconda container:
#
FROM continuumio/anaconda3:latest

#
# We'll need development tools:
#
RUN apt-get -y update
RUN apt-get -y install build-essential
RUN conda update -n base -c defaults conda

#
# Install NumPy, Cython, and Python 3.8 or better into
# a new environment and set that as the default in root's
# .bashrc:
#
RUN conda install numpy cython>=0.22

#
# Copy CKBIT source  into the container:
#
RUN mkdir /tmp/ckbit
COPY pystan-install.sh pystan-test.py setup.py /tmp/ckbit/
COPY ckbit /tmp/ckbit/ckbit
WORKDIR /tmp/ckbit

#
# Build and install the CVODES-enabled PyStan:
#
RUN chmod +x pystan-install.sh
RUN ./pystan-install.sh --build-root /tmp/ckbit --verbose
 
#
# Install CKBIT:
#
RUN python setup.py install

#
# Remove the CKBIT source:
#
#WORKDIR /
#RUN rm -rf /tmp/ckbit
