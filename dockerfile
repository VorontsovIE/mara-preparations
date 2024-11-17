# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory in the container
WORKDIR /work

# Copy the application code into the container
# COPY src/ /work
COPY src/MaraPreprocessingApp.py /work/src/
COPY config.json /work
COPY source_data/ /work
COPY tss_clusters.bed /work

# Install required system packages
RUN apt-get update && apt-get install -y \
    wget \
    default-jre \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Corrected ENTRYPOINT syntax
# ENTRYPOINT ["python", "/app/src/MaraPreprocessingApp.py"]