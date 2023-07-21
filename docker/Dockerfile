FROM python:3.9.17-slim-bookworm

MAINTAINER Ben Weisburd

# install java-8    (adds ~350Mb)
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-transport-https \
    ca-certificates \
    wget
RUN mkdir -p /etc/apt/keyrings \
    && wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc \
    && echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list;

RUN apt-get update \
    && apt-get install -y --no-install-recommends temurin-8-jdk

RUN apt-get update && apt-get install -y --no-install-recommends \
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
    libncursesw5-dev \
    libgdbm-dev \
    libc6-dev \
    zlib1g-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    libffi-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    procps \
    less \
    htop \
    vim

ENV SAMTOOLS_VERSION="1.9"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN apt-get install -y --no-install-recommends libncursesw5-dev
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

RUN samtools --version

# python3 packages
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install --upgrade tqdm pyfaidx pybedtools cyvcf2 peewee scipy numpy pandas matplotlib seaborn ipython

# install gcloud
ENV GCLOUD_SDK_VERISON="438.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# install sqlite3
RUN wget https://www.sqlite.org/2023/sqlite-autoconf-3420000.tar.gz \
    && tar xzf sqlite-autoconf-3420000.tar.gz \
    && cd sqlite-autoconf-3420000 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf sqlite-autoconf-3330000*

RUN python3 -m pip install hail crcmod

# install cloudfuse
RUN GCSFUSE_REPO=gcsfuse-bullseye \
	&& echo "deb https://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list \
	&& apt-get install --no-install-recommends --yes gnupg \
	&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
	&& apt-get update && apt-get install --no-install-recommends --yes fuse gcsfuse


# install all the GATK versions used in generating the gnomAD v4 callset
COPY --from=broadinstitute/gatk:4.0.10.1 /gatk/gatk-package-4.0.10.1-local.jar /gatk/gatk-v4.0.10.1.jar
COPY --from=broadinstitute/gatk:4.1.0.0 /gatk/gatk-package-4.1.0.0-local.jar /gatk/gatk-v4.1.0.0.jar
COPY --from=broadinstitute/gatk:4.1.4.1 /gatk/gatk-package-4.1.4.1-local.jar /gatk/gatk-v4.1.4.1.jar
COPY --from=broadinstitute/gatk:4.1.8.0 /gatk/gatk-package-4.1.8.0-local.jar /gatk/gatk-v4.1.8.0.jar

COPY --from=weisburd/base-other-tools:latest /usr/local/bin/bedtools /usr/local/bin/bedtools
RUN bedtools --version

ENTRYPOINT ["/bin/bash"]
