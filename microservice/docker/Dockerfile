FROM rocker/r-ver:4.3.2

ENV R_LIBS_USER=/usr/local/lib/R/site-library \
    R_LIBS_SITE=/usr/local/lib/R/site-library \
    PLUMBER_PORT=8002 \
    PLUMBER_HOST=0.0.0.0 \
    APP_ROOT=/app \
    TF_CPP_MIN_LOG_LEVEL=2

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    git \
    pkg-config \
    cmake \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libicu-dev \
    zlib1g-dev \
    libopenblas-dev \
    liblapack-dev \
    libgit2-dev \
    libsodium-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libcairo2-dev \
    curl \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Ensure R library directory exists and is used for installations
RUN mkdir -p /usr/local/lib/R/site-library && chmod 777 /usr/local/lib/R/site-library

# Prefer Posit Package Manager binaries when available (faster, more robust)
ENV RSPM=https://packagemanager.posit.co/cran/__linux__/bookworm/latest
RUN bash -lc "echo \"options(repos=c(RSPM='${RSPM}', CRAN='https://cran.r-project.org'))\" >> /usr/local/lib/R/etc/Rprofile.site && cat /usr/local/lib/R/etc/Rprofile.site"

# Install required R packages individually to handle dependencies better
RUN R -e "install.packages('jsonlite', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages('stringi', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages('plumber', lib='/usr/local/lib/R/site-library', dependencies=TRUE); if (!requireNamespace('plumber', quietly=TRUE, lib.loc='/usr/local/lib/R/site-library')) { message('plumber not available after install'); q(status=1) }"
# seqinr needs zlib (already installed) and imports ade4, segmented
RUN R -e "install.packages(c('ade4','segmented'), lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "if (!requireNamespace('seqinr', quietly=TRUE, lib.loc='/usr/local/lib/R/site-library')) { install.packages('seqinr', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE) } ; if (!requireNamespace('seqinr', quietly=TRUE, lib.loc='/usr/local/lib/R/site-library')) { q(status=1) }"
RUN R -e "install.packages('BiocManager', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/'); BiocManager::install(c('protr'), ask=FALSE, update=FALSE, lib='/usr/local/lib/R/site-library')"
RUN R -e "install.packages('randomForest', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages('ROCR', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"
## Pre-install caret heavy dependencies to avoid long chained builds
RUN R -e "install.packages(c('Rcpp','RcppEigen','Matrix','minqa','nloptr'), lib='/usr/local/lib/R/site-library', dependencies=TRUE)"
RUN R -e "install.packages(c('lme4','ggplot2','lattice','foreach','ModelMetrics','pROC','recipes','withr','generics','R6'), lib='/usr/local/lib/R/site-library', dependencies=TRUE)"
RUN R -e "options(Ncpus=max(1, parallel::detectCores()-1)); install.packages('caret', lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"
RUN R -e "if (!requireNamespace('caret', quietly=TRUE, lib.loc='/usr/local/lib/R/site-library')) { message('caret not available after install'); q(status=1) }"
RUN R -e "install.packages(c('reticulate','keras'), lib='/usr/local/lib/R/site-library', repos='https://cran.r-project.org/', dependencies=TRUE)"

# Python ML runtime for reticulate-backed Keras (TF 2.13 + h5py)
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel && \
    python3 -m pip install --no-cache-dir \
      tensorflow-aarch64==2.13.0 \
      h5py==3.10.0 \
      numpy==1.24.3

# Bind R reticulate to system python by default
ENV RETICULATE_PYTHON=/usr/bin/python3

# Ensure runtime library path includes our site-library
ENV R_LIBS=/usr/local/lib/R/site-library

# Diagnostics: show R library paths and what packages are installed where
RUN R -q -e "print(.libPaths()); cat('\nContents of /usr/local/lib/R/site-library:\n'); print(list.files('/usr/local/lib/R/site-library'))"

# Verify each package can be loaded
RUN R -e "cat('Testing package loading...\n'); library(jsonlite); cat('jsonlite OK\n')"
RUN R -e "library(plumber); cat('plumber OK\n')"
RUN R -e "library(seqinr); cat('seqinr OK\n')"
RUN R -e "library(protr); cat('protr OK\n')"
RUN R -e "library(randomForest); cat('randomForest OK\n')"
RUN R -e "library(ROCR); cat('ROCR OK\n')"
RUN R -e "if (!requireNamespace('caret', quietly=TRUE, lib.loc='/usr/local/lib/R/site-library')) { q(status=1) } ; library(caret); cat('caret OK\n')"
RUN R -e "library(reticulate); cat('reticulate OK -> using Python: ', reticulate::py_config()\$python, '\n')"
RUN R -e "library(keras); cat('keras OK (R) -> TF available: ', tryCatch(keras::is_keras_available(), error=function(e) FALSE), '\n')"

WORKDIR /app

# Copy project files
COPY . /app/

EXPOSE 8002

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD curl -f http://localhost:8002/health || exit 1

CMD ["bash", "microservice/scripts/start.sh"]

