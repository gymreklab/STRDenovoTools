FROM gymreklab/str-toolkit:latest

# Get necessary packages
RUN apt-get install -qqy \
  autotools-dev \
  automake \
  libgsl-dev

# Download, compile, and install CookieMonSTR
RUN git clone https://github.com/gymreklab/STRDenovoTools
WORKDIR STRDenovoTools
RUN autoreconf --install
RUN autoconf
RUN ./configure
RUN make
RUN make install
WORKDIR ..
