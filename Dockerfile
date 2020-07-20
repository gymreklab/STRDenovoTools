FROM ileena/strtools:latest

# Get necessary packages
RUN apt-get update

RUN apt-get install -qqy \
  autotools-dev \
  automake \
  libtool libtool-bin \
  libgsl-dev


#install latest cmake (needed for STRDenovoTools)
ADD https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh /cmake-3.7.2-Linux-x86_64.sh
RUN mkdir /opt/cmake
RUN sh /cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license
RUN ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
RUN cmake --version


# Install packages
RUN pip3 install pycrypto
RUN pip3 install pandas


# Download, compile, and install GangSTR v2.4.4 for chrX
RUN wget -O GangSTR-2.4.4.tar.gz https://github.com/gymreklab/GangSTR/releases/download/v2.4.4/GangSTR-2.4.4.tar.gz
RUN tar -xzvf GangSTR-2.4.4.tar.gz
WORKDIR GangSTR-2.4.4
RUN ./install-gangstr.sh
RUN ldconfig
WORKDIR ..



# Download, compile, and install Vcftools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
RUN tar -xzvf vcftools-0.1.16.tar.gz
WORKDIR vcftools-0.1.16
RUN ./autogen.sh && ./configure && make && make install
WORKDIR ..


# Download, compile, and install CookieMonSTR
RUN apt-get update && apt-get install -qqy cmake
RUN git clone https://github.com/ileenamitra/STRDenovoTools
WORKDIR STRDenovoTools
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make
RUN make install
WORKDIR ..
WORKDIR ..


# Download, compile, and install Datamash
RUN wget http://ftp.gnu.org/gnu/datamash/datamash-1.3.tar.gz
RUN tar -xzf datamash-1.3.tar.gz
WORKDIR datamash-1.3
RUN ./configure
RUN make
RUN make check
RUN make install
