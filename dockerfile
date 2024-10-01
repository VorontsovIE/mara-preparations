# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the application code into the container
COPY app/ /app/app/
COPY config.json /app/config.json
COPY data/ /app/data/
COPY MaraPreprocessingApp.py /app/MaraPreprocessingApp.py

# Install required system packages
RUN apt-get update && apt-get install -y     wget     default-jre     bedtools     && rm -rf /var/lib/apt/lists/*

# Set the entry point to the application
ENTRYPOINT [python, app/MaraPreprocessingApp.py]
