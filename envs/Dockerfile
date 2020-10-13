FROM continuumio/miniconda:4.7.12

WORKDIR /home/snakemake/

COPY environment.yaml ./


# mamba is a faster C++ re-implemenbtation of conda
RUN conda install -c conda-forge mamba --yes \
  && mamba create -c bioconda \
                  -c conda-forge \
                  -c r \
                  --name kallisto \
                  r-dplyr=1.0.2 \
                  fastp=0.19.5 \
                  kallisto=0.45.0 \
                  r-optparse=1.6.6 \
                  r-sleuth=0.30.0 \
                  snakemake=5.26.0 \
  && conda clean --all

RUN echo "source activate kallisto" > ~/.bashrc
ENV PATH /opt/conda/envs/kallisto/bin:$PATH

ENTRYPOINT ["snakemake"]