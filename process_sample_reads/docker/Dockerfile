FROM bitnami/minideb:stretch

MAINTAINER Ben Weisburd

# install basic linux utils like 'uptime' and 'free'
RUN install_packages procps

# install java-8    (adds ~350Mb)
RUN install_packages default-jdk

RUN apt-get update \
	&& apt-get dist-upgrade -y \
	&& apt-get install -y python-smbus libncursesw5-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev

RUN install_packages \
    ca-certificates \
    wget \
    curl \
    pkg-config \
    gcc \
    man-db \
    g++ \
    make \
    autoconf \
    unzip \
    bzip2 \
    zlib1g-dev \
    dpkg-dev \
    build-essential \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev

ENV PYTHON_VERSION="3.7.1"
RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
	&& tar xvf Python-${PYTHON_VERSION}.tar.xz \
	&& rm Python-${PYTHON_VERSION}.tar.xz \
	&& cd Python-${PYTHON_VERSION} \
	&& ./configure --enable-optimizations \
	&& make install \
	&& make clean

ENV SAMTOOLS_VERSION="1.9"

RUN install_packages wget libbz2-dev zlib1g-dev liblzma-dev libcurl4-openssl-dev build-essential
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN install_packages libncursesw5-dev
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN samtools --version

COPY --from=broadinstitute/gatk:4.1.4.1 /gatk/gatk-package-4.1.4.1-local.jar /gatk/gatk.jar
RUN java -jar /gatk/gatk.jar --version

COPY --from=us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135 /usr/gitc/GATK35.jar /gatk/GATK35.jar
RUN java -jar /gatk/GATK35.jar --version

COPY --from=weisburd/base-other-tools:latest /usr/local/bin/bedtools /usr/local/bin/bedtools
RUN bedtools --version

# python3 packages
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install --upgrade tqdm pyvcf pyfaidx pybedtools cyvcf2 peewee scipy numpy pandas matplotlib seaborn ipython

# install gcloud
ENV GCLOUD_SDK_VERISON="297.0.1"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH


# install hail
RUN python3 -m pip install hail

# nice-to-have linux utils
RUN install_packages less htop vim emacs
