FROM ubuntu:disco
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --install-recommends apt-utils software-properties-common dirmngr apt-transport-https build-essential 
RUN apt-get install -y python3.7 python3-pip python3-setuptools python3-dev
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/'
RUN apt-get update && apt-get install -y --install-recommends \
default-jdk \
nano \
libblas-dev \
liblapack-dev \
gfortran \
libgmp3-dev \
libcurl4-openssl-dev \
libxml2-dev \
r-base \
r-base-dev \
r-cran-matrix         

RUN su - -c "R CMD javareconf"

WORKDIR /app
COPY R/install.packages.R /app/R/install.packages.R

RUN echo ### Install R packages ###
RUN Rscript R/install.packages.R


COPY . /app
