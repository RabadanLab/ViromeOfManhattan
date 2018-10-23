############################################################
# Dockerfile to run Pandora
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER docliver

################## INSTALL STUFF ######################
# Install programs

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
RUN apt-get install -y cmake
RUN apt-get install -y s3fs

# Install python packages
RUN pip install awscli
RUN pip install numpy scipy matplotlib pandas sympy nose biopython

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

# RUN https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
# RUN tar -zxvf jellyfish-2.2.10.tar.gz
# WORKDIR /home/jellyfish-2.2.10
# RUN ./configure && make && make install
# WORKDIR /home

# RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.11.3/salmon-0.11.3-linux_x86_64.tar.gz
# RUN tar -zxvf salmon-0.11.3-linux_x86_64.tar.gz

# RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.4.tar.gz
# RUN tar -zxvf Trinity-v2.8.4.tar.gz
# WORKDIR /home/trinityrnaseq-Trinity-v2.8.4
# RUN make
# WORKDIR /home

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz
RUN tar -zxvf v2.1.1.tar.gz
WORKDIR /home/trinityrnaseq-2.1.1
RUN rm Chrysalis/Makefile
COPY Chrysalis/Makefile Chrysalis/Makefile
RUN make
WORKDIR /home

RUN wget https://sourceforge.net/projects/subread/files/subread-1.6.2/subread-1.6.2-Linux-x86_64.tar.gz
RUN tar -zxvf subread-1.6.2-Linux-x86_64.tar.gz
RUN git clone https://github.com/RabadanLab/Pandora.git
# BEWARE: docker build relies on cache rather than grabbing repo updates
# WORKDIR /home/Pandora
# RUN git pull
# WORKDIR /home
# RUN rm *.gz *.tar *.zip'

RUN mkdir -p data results

# ADD COMMANDS TO PATH
ENV PATH="${PATH}:/home/STAR-2.6.0a/bin/Linux_x86_64"
ENV PATH="${PATH}:/home/bowtie2-2.3.4.3-linux-x86_64"
ENV PATH="${PATH}:/home/ncbi-blast-2.7.1+/bin"
ENV PATH="${PATH}:/home/subread-1.6.2-Linux-x86_64/bin"
# ENV PATH="${PATH}:/home/trinityrnaseq-Trinity-v2.8.4"
ENV PATH="${PATH}:/home/trinityrnaseq-2.1.1"
# ENV PATH="${PATH}:/home/salmon-0.11.3-linux_x86_64/bin"

##################### INSTALL STUFF #####################

# SET UP DIRECTORIES, COPY FILES
RUN mkdir -p data ref results /opt/software /root/.aws
RUN ln -s /home/Pandora /opt/software/Pandora

# COPY DATA
COPY data/mate1.fq.gz /home/data/
COPY data/mate2.fq.gz /home/data/

# RUN PANDORA
WORKDIR /home/results
# CMD echo "hello world"
CMD python ../Pandora/pandora.py scan -id 1 --verbose --gzip -c ../Pandora/pandora.config.docker.txt -r1 ../data/mate1.fq.gz -r2 ../data/mate2.fq.gz --noclean --verbose --noSGE
