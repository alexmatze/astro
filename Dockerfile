#
# base
#
FROM  ubuntu:16.04

#
# common-environment
#
ENV USER lofar
ENV INSTALLDIR /home/${USER}/opt

#
# environment
#
ENV DEBIAN_FRONTEND noninteractive
ENV PYTHON_VERSION 2.7

#
# versions
#
ENV CFITSIO_VERSION 3390
ENV WCSLIB_VERSION 5.15
ENV LOG4CPLUS_VERSION 1.1.x
ENV CASACORE_VERSION v2.2.0
ENV CASAREST_VERSION v1.4.1
ENV PYTHON_CASACORE_VERSION v2.1.2
ENV AOFLAGGER_VERSION v2.8.0
ENV LOFAR_VERSION 2_19_5

#
# set-uid
#
ENV UID 1000

#
# build environment
#
ENV J 2

#
# base
#
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get -y install sudo
RUN apt-get -y install apt-utils
RUN apt-get -y install git subversion wget
RUN apt-get -y install automake autotools-dev cmake make python-setuptools
RUN apt-get -y install  g++ gcc gfortran
RUN apt-get -y install libblas-dev libfftw3-dev python-dev liblapack-dev libpng-dev libxml2-dev python-numpy libreadline-dev libncurses-dev python-scipy liblog4cplus-dev
RUN apt-get -y install libboost-dev libboost-python-dev libboost-thread-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev libboost-signals-dev
RUN apt-get -y install bison bzip2 flex python-xmlrunner python-pip gettext doxygen libgsl-dev libhdf5-dev libboost-test-dev
RUN pip install --upgrade pip
RUN pip install pyfits pywcs python-monetdb unittest2 matplotlib astropy aplpy
RUN cd /usr/lib/python2.7/dist-packages/numpy/core && ln -s `ls multiarray.*.so | head`  multiarray.so

#
# Additional stuff
#
RUN apt-get -y install ssh
RUN apt-get -y install vim
RUN apt-get -y install net-tools
RUN apt-get -y install htop
RUN apt-get -y install screen

RUN service ssh restart
RUN echo "defshell -bash" > ~/.screenrc





#
# setup-account
#
RUN getent group sudo &>/dev/null || groupadd sudo
RUN echo "useradd -m ${USERADD_FLAGS} ${USER}"
RUN useradd -m -u ${UID} ${USER}
RUN usermod -a -G sudo ${USER}
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN sed -i 's/requiretty/!requiretty/g' /etc/sudoers


USER ${USER}

#
# install-cfitsio
#
RUN mkdir -p ${INSTALLDIR}/cfitsio/build
RUN cd ${INSTALLDIR}/cfitsio && wget --retry-connrefused ftp://anonymous@heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio${CFITSIO_VERSION}.tar.gz
RUN cd ${INSTALLDIR}/cfitsio && tar xf cfitsio${CFITSIO_VERSION}.tar.gz
RUN cd ${INSTALLDIR}/cfitsio/build && cmake -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/cfitsio/ ../cfitsio
RUN cd ${INSTALLDIR}/cfitsio/build && make -j ${J}
RUN cd ${INSTALLDIR}/cfitsio/build && make install

#
# install-wcslib
#
RUN mkdir ${INSTALLDIR}/wcslib
RUN if [ "${WCSLIB_VERSION}" = "latest" ]; then cd ${INSTALLDIR}/wcslib && wget --retry-connrefused ftp://anonymous@ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2 -O wcslib-latest.tar.bz2; fi
RUN if [ "${WCSLIB_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/wcslib && wget --retry-connrefused ftp://anonymous@ftp.atnf.csiro.au/pub/software/wcslib/wcslib-${WCSLIB_VERSION}.tar.bz2; fi
RUN cd ${INSTALLDIR}/wcslib && tar xf wcslib-*.tar.bz2
RUN cd ${INSTALLDIR}/wcslib/wcslib* && ./configure --prefix=${INSTALLDIR}/wcslib --with-cfitsiolib=${INSTALLDIR}/cfitsio/lib/ --with-cfitsioinc=${INSTALLDIR}/cfitsio/include/ --without-pgplot
RUN cd ${INSTALLDIR}/wcslib/wcslib* && make
RUN cd ${INSTALLDIR}/wcslib/wcslib* && make install

#
# install-casacore
#
RUN mkdir -p ${INSTALLDIR}/casacore/build
RUN mkdir -p ${INSTALLDIR}/casacore/data
RUN cd ${INSTALLDIR}/casacore && git clone https://github.com/casacore/casacore.git src
RUN if [ "${CASACORE_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/casacore/src && git checkout tags/${CASACORE_VERSION}; fi
RUN cd ${INSTALLDIR}/casacore/data && wget --retry-connrefused ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar
RUN cd ${INSTALLDIR}/casacore/data && tar xf WSRT_Measures.ztar
RUN cd ${INSTALLDIR}/casacore/build && cmake -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/casacore/ -DDATA_DIR=${INSTALLDIR}/casacore/data -DWCSLIB_ROOT_DIR=/${INSTALLDIR}/wcslib/ -DCFITSIO_ROOT_DIR=${INSTALLDIR}/cfitsio/ -DBUILD_PYTHON=True -DUSE_OPENMP=True -DUSE_FFTW3=TRUE ../src/
RUN cd ${INSTALLDIR}/casacore/build && make -j ${J}
RUN cd ${INSTALLDIR}/casacore/build && make install

#
# install-casarest
#
RUN mkdir -p ${INSTALLDIR}/casarest/build
RUN cd ${INSTALLDIR}/casarest && git clone https://github.com/casacore/casarest.git src
RUN if [ "${CASAREST_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/casarest/src && git checkout tags/${CASAREST_VERSION}; fi
RUN cd ${INSTALLDIR}/casarest/build && cmake -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/casarest -DCASACORE_ROOT_DIR=${INSTALLDIR}/casacore -DWCSLIB_ROOT_DIR=${INSTALLDIR}/wcslib -DCFITSIO_ROOT_DIR=${INSTALLDIR}/cfitsio ../src/
RUN cd ${INSTALLDIR}/casarest/build && make -j ${J}
RUN cd ${INSTALLDIR}/casarest/build && make install

#
# install-python-casacore
#
RUN mkdir ${INSTALLDIR}/python-casacore
RUN cd ${INSTALLDIR}/python-casacore && git clone https://github.com/casacore/python-casacore
RUN if [ "$PYTHON_CASACORE_VERSION" != "latest" ]; then cd ${INSTALLDIR}/python-casacore/python-casacore && git checkout tags/${PYTHON_CASACORE_VERSION}; fi
RUN cd ${INSTALLDIR}/python-casacore/python-casacore && ./setup.py build_ext -I${INSTALLDIR}/wcslib/include:${INSTALLDIR}/casacore/include/:${INSTALLDIR}/cfitsio/include -L${INSTALLDIR}/wcslib/lib:${INSTALLDIR}/casacore/lib/:${INSTALLDIR}/cfitsio/lib/ -R${INSTALLDIR}/wcslib/lib:${INSTALLDIR}/casacore/lib/:${INSTALLDIR}/cfitsio/lib/
RUN mkdir -p ${INSTALLDIR}/python-casacore/lib/python${PYTHON_VERSION}/site-packages/
RUN mkdir -p ${INSTALLDIR}/python-casacore/lib64/python${PYTHON_VERSION}/site-packages/
RUN export PYTHONPATH=${INSTALLDIR}/python-casacore/lib/python${PYTHON_VERSION}/site-packages:${INSTALLDIR}/python-casacore/lib64/python${PYTHON_VERSION}/site-packages:$PYTHONPATH && cd ${INSTALLDIR}/python-casacore/python-casacore && ./setup.py develop --prefix=${INSTALLDIR}/python-casacore/


#
# install-aoflagger
#
RUN mkdir -p ${INSTALLDIR}/aoflagger/build
RUN cd ${INSTALLDIR}/aoflagger && git clone git://git.code.sf.net/p/aoflagger/code aoflagger
RUN cd ${INSTALLDIR}/aoflagger/aoflagger && git checkout tags/${AOFLAGGER_VERSION}
RUN cd ${INSTALLDIR}/aoflagger/build && cmake -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/aoflagger/ -DCASACORE_ROOT_DIR=${INSTALLDIR}/casacore -DCFITSIO_ROOT_DIR=${INSTALLDIR}/cfitsio -DBUILD_SHARED_LIBS=ON ../aoflagger
RUN cd ${INSTALLDIR}/aoflagger/build && make -j ${J}
RUN cd ${INSTALLDIR}/aoflagger/build && make install

#
# install-log4cplus
#
RUN echo "log4cplus is installed from the repositories"

#
# install-lofar
#
RUN mkdir -p ${INSTALLDIR}/lofar/build/gnu_opt
RUN if [ "${LOFAR_VERSION}" = "latest" ]; then cd ${INSTALLDIR}/lofar && svn --non-interactive -q co https://svn.astron.nl/LOFAR/trunk src; fi
RUN if [ "${LOFAR_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/lofar && svn --non-interactive -q co https://svn.astron.nl/LOFAR/tags/LOFAR-Release-${LOFAR_VERSION} src; fi
RUN cd ${INSTALLDIR}/lofar/build/gnu_opt && cmake -DBUILD_PACKAGES=Offline -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/lofar/ -DWCSLIB_ROOT_DIR=${INSTALLDIR}/wcslib/ -DCFITSIO_ROOT_DIR=${INSTALLDIR}/cfitsio/ -DCASAREST_ROOT_DIR=${INSTALLDIR}/casarest/ -DCASACORE_ROOT_DIR=${INSTALLDIR}/casacore/ -DAOFLAGGER_ROOT_DIR=${INSTALLDIR}/aoflagger/ -DLOG4CPLUS_ROOT_DIR=${INSTALLDIR}/log4cplus/ -DUSE_OPENMP=True ${INSTALLDIR}/lofar/src/
RUN cd ${INSTALLDIR}/lofar/build/gnu_opt && make -j ${J}
RUN cd ${INSTALLDIR}/lofar/build/gnu_opt && make install


#
# Fix Python .egg folders
#



#
# WSClean install
# 
RUN sudo mkdir -p ${INSTALLDIR}/wsclean
RUN sudo chmod -R 777 ${INSTALLDIR}/wsclean
RUN cd ${INSTALLDIR}/wsclean && wget https://sourceforge.net/projects/wsclean/files/wsclean-2.4/wsclean-2.4.tar.bz2
RUN sudo tar xf ${INSTALLDIR}/wsclean/wsclean-2.4.tar.bz2
RUN sudo sh -c 'echo ls'
RUN sudo rm ${INSTALLDIR}/wsclean/wsclean-2.4.tar.bz2
RUN sudo chmod -R 777 ${INSTALLDIR}/wsclean/wsclean-2.4
RUN cd ${INSTALLDIR}/wsclean/wsclean-2.4
RUN mkdir build && cd build
RUN cmake -DCMAKE_PREFIX_PATH="/home/lofar/opt/lofar/;/home/lofar/opt/casacore/;/home/lofar/opt/cfitsio;/home/lofar/opt/lofar/lib/python2.7/site-packages/losoto/" -DBUILD_SHARED_LIBS=TRUE ../
RUN make -j ${J}
RUN make install



#
# init-lofar
#
RUN sudo sh -c 'echo source \${INSTALLDIR}/lofar/lofarinit.sh  >> /usr/bin/init-lofar.sh'
RUN sudo sh -c 'echo export PYTHONPATH=\${PYTHONPATH:+:\${PYTHONPATH}}:\${INSTALLDIR}/python-casacore/lib/python2.7/site-packages/  >> /usr/bin/init-lofar.sh'
RUN sudo sh -c 'echo export PATH=\${PATH:+:\$PATH}:\${INSTALLDIR}/casacore/bin  >> /usr/bin/init-lofar.sh'
RUN sudo sh -c 'echo export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}:\${INSTALLDIR}/casacore/lib  >> /usr/bin/init-lofar.sh'
RUN sudo sh -c "echo source /usr/bin/init-lofar.sh >> /usr/bin/init.sh"





#
# entrypoint
#
ENTRYPOINT /bin/bash --init-file /usr/bin/init.sh

