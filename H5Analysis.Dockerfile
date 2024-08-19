FROM jupyter/datascience-notebook@sha256:9504f4f4ab7e89b49d61d7be2e9ff8c57870de2050aa4360f55b2e59193f7486 as h5analysis

WORKDIR /home

COPY h5analysis_requirements.txt .
RUN pip install -r h5analysis_requirements.txt

EXPOSE 8888
