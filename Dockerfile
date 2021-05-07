#
# This is the definition of the CKBIT Docker container.  Derives from
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
# Install NumPy and Cython into the default environment:
#
RUN conda install numpy cython>=0.22

#
# Copy CKBIT source into the container:
#
RUN mkdir /tmp/ckbit
COPY pystan-install.sh pystan-test.py setup.py /tmp/ckbit/
COPY ckbit /tmp/ckbit/ckbit
WORKDIR /tmp/ckbit

#
# Build and install the CVODES-enabled PyStan: and install CKBIT: and remove build folders
# This reduces image size by reducing the last layers. 
#
RUN chmod +x pystan-install.sh
RUN ./pystan-install.sh --build-root /tmp/ckbit --verbose && python setup.py install && rm -rf /tmp/ckbit
WORKDIR /

