FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y alien rpm && rm -rf /var/lib/apt/lists/*
COPY bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm /tmp/bcl2fastq.rpm
RUN alien -i /tmp/bcl2fastq.rpm && rm /tmp/bcl2fastq.rpm
RUN bcl2fastq --version
CMD ["bash"]