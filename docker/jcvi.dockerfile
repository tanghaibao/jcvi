FROM continuumio/miniconda3
ARG DEBIAN_FRONTEND=noninteractive
LABEL maintainer="tanghaibao@gmail.com"

RUN apt-get update
RUN apt-get install -y gcc git build-essential
RUN apt-get install -y python3-dev libxml2-dev libxslt-dev
RUN apt-get install -y libncurses-dev libcurl4-openssl-dev zlib1g-dev
RUN apt-get install -y vcftools python3-pip
RUN apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
RUN apt-get install -y wget autoconf libssl-dev

RUN pip install boto3 awscli
RUN pip install pyfaidx pyliftover
RUN pip install cython
RUN pip install pandas
RUN pip install scipy

# GRABSEEDS dependencies
# <https://github.com/tanghaibao/jcvi/wiki/GRABSEEDS%3A-How-to-install>
RUN apt-get install -y libxft-dev libfreetype6 libfreetype6-dev
RUN apt-get install -y libmagickwand-dev
RUN apt-get install -y texlive texlive-latex-extra texlive-latex-recommended cm-super
RUN apt-get install -y dvipng

RUN pip install matplotlib scikit-image pypdf2 wand Pillow

# https://github.com/tanghaibao/jcvi/issues/307
RUN apt install curl

# Install jcvi
# https://github.com/tanghaibao/jcvi/issues/509
RUN pip install jcvi
WORKDIR /
