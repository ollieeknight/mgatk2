FROM python:3.10-slim

WORKDIR /opt/mgatk2
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir . && mgatk2 paired --help
