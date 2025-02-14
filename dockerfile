FROM python:3.9-slim

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /work

# COPY src/ /work
COPY src/MaraPreprocessingApp.py /work/src/
COPY config.json /work
COPY source_data/ /work/source_data/
COPY tss_clusters.bed /work

RUN apt-get update && apt-get install -y \
    wget \
    default-jre \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# ENTRYPOINT ["python", "/app/src/MaraPreprocessingApp.py"]
