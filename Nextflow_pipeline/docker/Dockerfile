FROM mambaorg/micromamba:2.0.8

USER $MAMBA_USER

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    r-devtools \
    r-r.oo=1.25.0 \
    r-gtools=3.9.4 \
    r-stringr=1.5.0 \
    r-purrr=1.0.1 \
    r-colorspace=2.1_0 \
    r-r.methodss3=1.8.2 \
    r-seuratobject=5.0.2 \
    r-tidyr=1.3.0 \
    r-r.utils=2.12.2 \
    r-dplyr=1.1.4 \
    r-ggplot2=3.5.2 \
    r-seurat=4.3.0.1 \
    r-missMDA=1.18 \
    r-dismo \
    r-dbscan \
    bioconductor-biocparallel=1.34.2 \
    bioconductor-complexheatmap && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN which R

RUN R -e 'devtools::install_github("crazyhottommy/scclusteval", dependencies=FALSE)'

ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
