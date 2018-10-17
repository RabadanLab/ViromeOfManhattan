############################################################
# Dockerfile to run Pandora
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER docliver

################## BEGIN INSTALLATION ######################
# Install 

# Usage: WORKDIR /path
WORKDIR /home

# Update the repository sources list
RUN apt-get update

# Install ubuntu packages
RUN apt-get install -y automake autotools-dev fuse g++ git libfuse-dev libssl-dev libxml2-dev make pkg-config
RUN apt-get install -y wget
RUN apt-get install -y gcc
RUN apt-get install -y libgomp1
RUN apt-get install -y perl
RUN apt-get install -y python
RUN apt-get install -y libidn11
RUN apt-get install -y curl
RUN apt-get install -y perl 
RUN apt-get install -y python 
RUN apt-get install -y libidn11 
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libcurl4-gnutls-dev
RUN apt-get install -y build-essential
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y libpng-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev
RUN apt-get install -y python-pip
RUN apt-get install -y unzip
RUN apt-get install -y rsync
RUN apt-get install -y s3fs

# Install python packages
RUN pip install awscli
RUN pip install numpy scipy matplotlib pandas sympy nose biopython

###

# Install (and build) bioinformatics packages
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
RUN tar -zxvf ncbi-blast-2.7.1+-x64-linux.tar.gz
RUN wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
RUN tar -zxvf 2.6.0a.tar.gz
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.4.1/samtools-1.4.1.tar.bz2
RUN bunzip2 samtools-1.4.1.tar.bz2
RUN tar -xvf samtools-1.4.1.tar
WORKDIR /home/samtools-1.4.1
RUN ./configure && make && make install
WORKDIR /home
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.4.1/htslib-1.4.1.tar.bz2
RUN bunzip2 htslib-1.4.1.tar.bz2
RUN tar -xvf htslib-1.4.1.tar
WORKDIR /home/htslib-1.4.1
RUN ./configure && make && make install
WORKDIR /home
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip
RUN unzip bowtie2-2.3.4.3-linux-x86_64.zip
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.4.tar.gz
RUN tar -zxvf Trinity-v2.8.4.tar.gz
RUN wget https://sourceforge.net/projects/subread/files/subread-1.6.2/subread-1.6.2-Linux-x86_64.tar.gz
RUN tar -zxvf subread-1.6.2-Linux-x86_64.tar.gz
RUN git clone https://github.com/RabadanLab/Pandora.git
# RUN rm *.gz *.tar *.zip

##################### INSTALLATION END #####################

# ADD TO PATH
ENV PATH="${PATH}:/home/STAR-2.6.0a/bin/Linux_x86_64"
ENV PATH="${PATH}:/home/bowtie2-2.3.4.3-linux-x86_64"
ENV PATH="${PATH}:/home/ncbi-blast-2.7.1+/bin"
ENV PATH="${PATH}:/home/subread-1.6.2-Linux-x86_64/bin"
ENV PATH="${PATH}:/home/trinityrnaseq-Trinity-v2.8.4"

# MISC (make links)
RUN mkdir -p /opt/software
WORKDIR /opt/software
RUN ln -s /home/Pandora
WORKDIR /home

# COMMAND
CMD echo "Hello docker!"
