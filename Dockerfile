FROM ubuntu:focal
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/anaconda/bin:${PATH}"
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

#RUN sed -i -e 's|disco|eoan|g' /etc/apt/sources.list
RUN apt-get update && apt-get install -y --install-recommends apt-utils software-properties-common dirmngr gnupg apt-transport-https ca-certificates software-properties-common build-essential 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

#RUN apt-get update && apt-get install -y \
RUN apt-get update && apt-get install -y --install-recommends \
locales \
default-jdk \
nano \
vim \
wget \
zip \
unzip \
curl \
python2 \
r-base \
r-base-dev \
libcurl4-openssl-dev \
libssl-dev \
libssl1.1 \
libncurses5-dev

RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python2 get-pip.py

RUN echo "LC_ALL=en_US.UTF-8" >> /etc/environment
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen
RUN echo "LANG=en_US.UTF-8" > /etc/locale.conf
RUN locale-gen en_US.UTF-8

RUN su - -c "R CMD javareconf"

RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
RUN bash Anaconda3-2020.11-Linux-x86_64.sh -b -p /opt/anaconda
RUN rm Anaconda3-2020.11-Linux-x86_64.sh 

WORKDIR /CytoMeth
COPY . /CytoMeth

RUN yes | ./install.sh

COPY . /CytoMeth
