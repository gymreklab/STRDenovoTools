FROM gymreklab/str-toolkit-2.4:latest

# Get necessary packages
RUN apt-get install -qqy \
  autotools-dev \
  automake \
  libtool libtool-bin \
  libgsl-dev

# Download, compile, and install CookieMonSTR
RUN git clone https://github.com/gymreklab/STRDenovoTools
WORKDIR STRDenovoTools
RUN ./reconf
RUN ./configure
RUN make
RUN make install
WORKDIR ..
