FROM ubuntu:latest
MAINTAINER tanghaibao@gmail.com

RUN apt-get update
RUN apt-get install -y gcc git build-essential
RUN apt-get install -y python-dev libxml2-dev libxslt-dev
RUN apt-get install -y libncurses-dev libcurl4-openssl-dev zlib1g-dev
RUN apt-get install -y vcftools python-pip
RUN apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
RUN apt-get install -y wget autoconf libssl-dev

RUN pip install boto3 awscli
RUN pip install pyfaidx pyliftover pyvcf
RUN pip install cython
RUN pip install pandas
RUN pip install scipy

# GRABSEEDS dependencies
# <https://github.com/tanghaibao/jcvi/wiki/GRABSEEDS%3A-How-to-install>
RUN apt-get install -y libxft-dev libfreetype6 libfreetype6-dev
RUN apt-get install -y libmagickwand-dev
RUN apt-get install -y texlive texlive-latex-extra texlive-latex-recommended
RUN apt-get install -y dvipng

RUN pip install matplotlib scikit-image pypdf2 wand Pillow

# Install jcvi
RUN pip install git+git://github.com/tanghaibao/jcvi.git
WORKDIR /
