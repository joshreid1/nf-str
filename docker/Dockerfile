FROM mambaorg/micromamba:latest

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    libssl-dev \
    git \
    wget \
    libcurl4 \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libgomp1 \
    libgcc-s1 \
    libstdc++6 \
    procps \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /app
# Copy files with proper ownership
COPY --chown=$MAMBA_USER:$MAMBA_USER environments/ /app/environments/
RUN mkdir -p /usr/local/bin

# Install source packages via GitHub
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz \
    && tar -xvzf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz \
    && rm ExpansionHunter-v5.0.0-linux_x86_64.tar.gz \
    && mv ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter /usr/local/bin/ 

RUN wget https://github.com/G2Lab/scattr/releases/download/0.2.2/scattr-0.2.2-linux-x86_64.tar.gz \
    && tar -xvzf scattr-0.2.2-linux-x86_64.tar.gz \
    && rm scattr-0.2.2-linux-x86_64.tar.gz \
    && mv scattr /usr/local/bin/

RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v4.0.0/trgt-v4.0.0-x86_64-unknown-linux-gnu.tar.gz \
    && tar -xvzf trgt-v4.0.0-x86_64-unknown-linux-gnu.tar.gz \
    && rm trgt-v4.0.0-x86_64-unknown-linux-gnu.tar.gz \
    && mv trgt-v4.0.0-x86_64-unknown-linux-gnu/trgt /usr/local/bin/

RUN chmod +x /usr/local/bin/ExpansionHunter /usr/local/bin/scattr /usr/local/bin/trgt

RUN wget https://github.com/bahlolab/EH5_kmerfilter/archive/refs/tags/v1.0.0.tar.gz  \
    && tar -xvzf v1.0.0.tar.gz \
    && rm v1.0.0.tar.gz \
    && mv EH5_kmerfilter-1.0.0 /app/EH5_kmerfilter \
    && chmod +x /app/EH5_kmerfilter/kmer_filter.py \
    && ln -s /app/EH5_kmerfilter/kmer_filter.py /usr/local/bin/kmer_filter.py

USER $MAMBA_USER

RUN micromamba create -n eh5_kmerfilter -f /app/EH5_kmerfilter/environment.yml && \
    micromamba clean --all --yes 

# base dependencies for pipeline
RUN micromamba create -f /app/environments/base.yml -y && \
    micromamba clean --all --yes

# straglr
RUN micromamba create -f /app/environments/straglr.yml -y && \
    micromamba clean --all --yes

# ATaRVa
RUN micromamba create -f /app/environments/atarva.yml -y && \
    micromamba clean --all --yes

# longTR
RUN micromamba create -n longTR -c bioconda -c conda-forge longtr


ENV PATH="/usr/local/bin:${PATH}"
ENV ENV_NAME=base
