FROM ubuntu:disco
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --install-recommends apt-utils software-properties-common dirmngr apt-transport-https build-essential 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/'

#RUN apt-get update && apt-get install -y --install-recommends \
#RUN apt-get update && apt-get install -y \
default-jdk \
nano \
wget \
zip \
unzip \
python2\
r-base \
r-base-dev

RUN su - -c "R CMD javareconf"

RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
RUN bash Anaconda3-2019.03-Linux-x86_64.sh -b -p /opt/anaconda
RUN rm Anaconda3-2019.03-Linux-x86_64.sh 
RUN echo 'export PATH="/opt/anaconda/bin:$PATH"' >> ~/.bashrc
RUN source ~/.bashrc

WORKDIR /app
COPY . /app

RUN yes | ./install.sh

COPY . /app

