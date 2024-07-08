FROM continuumio/miniconda3:23.10.0-1

# Install dependencies
COPY ./environment.yml /tmp/environment.yml
RUN conda env update -n base --solver libmamba --file /tmp/environment.yml && conda clean -afy

WORKDIR /task