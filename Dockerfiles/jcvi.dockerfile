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

# Install jcvi
RUN pip install git+git://github.com/tanghaibao/jcvi.git
WORKDIR /
