FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        wget \
        git \
        unzip \
        procps \
        build-essential \
        openjdk-21-jre-headless \
        samtools \
        bcftools \
        tabix \
        bwa \
        minimap2 \
        fastqc \
        fastp \
        && \
    rm -rf /var/lib/apt/lists/*

COPY docker/mosdepth /usr/local/bin/mosdepth

RUN chmod +x /usr/local/bin/mosdepth && \
    mosdepth --version

COPY snpEff_latest_core.zip /tmp/snpEff_latest_core.zip

RUN mkdir -p /opt/snpEff && \
    unzip -q /tmp/snpEff_latest_core.zip -d /opt/snpEff && \
    rm /tmp/snpEff_latest_core.zip

ENV SNPEFF_HOME=/opt/snpEff/snpEff

RUN printf '#!/bin/sh\nexec java -jar /opt/snpEff/snpEff/snpEff.jar "$@"\n' \
    > /usr/local/bin/snpEff && \
    chmod +x /usr/local/bin/snpEff

RUN java -version && \
    snpEff -version

WORKDIR /opt/AmpSeeker

COPY requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir \
        jupyter \
        jupyterlab \
        "jupyter-book<2" \
        -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt


COPY fastplong /usr/local/bin/fastplong
RUN chmod +x /usr/local/bin/fastplong && fastplong --version

# Install kaleido and plotly (kaleido v1.x — no scope API)
RUN pip install --no-cache-dir "plotly>=5" "kaleido>=1.0"

# Install MultiQC via pip (>= 1.25 fixes missing JS assets in reports)
RUN pip install --no-cache-dir "multiqc>=1.25"

COPY . /opt/AmpSeeker

RUN python -m ipykernel install \
        --prefix=/usr/local \
        --name=AmpSeq_python \
        --display-name="AmpSeq_python"

RUN mkdir -p /opt/AmpSeeker/config \
             /opt/AmpSeeker/results \
             /opt/AmpSeeker/logs \
             /opt/AmpSeeker/resources \
             /opt/AmpSeeker/docs

CMD ["bash"]
