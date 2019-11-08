FROM ileena/strtools:latest

# Get necessary packages
RUN apt-get install -qqy \
  autotools-dev \
  automake \
  libtool libtool-bin \
  libgsl-dev


# Install packages
RUN pip3 install pycrypto
RUN pip3 install pandas==0.22.0


# Download, compile, and install GangSTR
RUN wget -O GangSTR-2.4.2.tar.gz https://github.com/gymreklab/GangSTR/releases/download/v2.4.2/GangSTR-2.4.2.tar.gz
RUN tar -xzvf GangSTR-2.4.2.tar.gz
WORKDIR GangSTR-2.4.2
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
RUN git clone https://github.com/gymreklab/STRDenovoTools
WORKDIR STRDenovoTools
RUN ./reconf
RUN ./configure
RUN make
RUN make install
WORKDIR ..

# Download, compile, and install Datamash
RUN wget http://ftp.gnu.org/gnu/datamash/datamash-1.3.tar.gz
RUN tar -xzf datamash-1.3.tar.gz
WORKDIR datamash-1.3
RUN ./configure
RUN make
RUN make check
RUN make install
