FROM bitnami/minideb:stretch

MAINTAINER Ben Weisburd

# install basic linux utils like 'uptime' and 'free'
RUN install_packages procps

# install java-8    (adds ~350Mb)
RUN install_packages default-jdk

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
