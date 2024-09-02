# LGM_H5N1
### General Dockerfile for the lgm_H5N1 project.
### Mostly based on jupyter datascience image.
### Uses a multi-stage build to add Jupyter
### capabilities to the nextstrain docker image.

# BUILD STAGE
FROM mambaorg/micromamba:1.5.7 AS build

## Install requirements from env.yml
USER root
COPY env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

# FINAL STAGE
FROM nextstrain/base as lgm_h5n1

## Copy the items installed by mamba during the build stage
COPY --from=build /opt/conda /opt/conda

## Parameters
ENV PATH="/opt/conda/bin:$PATH"
WORKDIR /home
EXPOSE 8888

## Run commands
### Launch Jupyter notebook server
CMD ["jupyter", "notebook", "--port=888", "--no-browser", "--allow-root"]
