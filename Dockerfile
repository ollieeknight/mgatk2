FROM python:3.10-slim

COPY . /opt/mgatk2
RUN pip install --no-cache-dir /opt/mgatk2
