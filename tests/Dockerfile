FROM r-base:4.3.2

# Install system dependencies needed by R packages (nloptr, car, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    cmake \
    libnlopt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R dependencies one-by-one to get clear error messages
RUN R -e "install.packages('reshape2', repos='https://cloud.r-project.org')" \
 && R -e "install.packages('nlme', repos='https://cloud.r-project.org')" \
 && R -e "install.packages('nloptr', repos='https://cloud.r-project.org')" \
 && R -e "install.packages('lme4', repos='https://cloud.r-project.org')" \
 && R -e "install.packages('car', repos='https://cloud.r-project.org')" \
 && R -e "install.packages('pbkrtest', repos='https://cloud.r-project.org')"

# Verify all deps are loadable
RUN R -e "library(reshape2); library(nlme); library(lme4); library(car); library(pbkrtest); cat('All dependencies OK\n')"

# Install msTrawler from source
COPY . /msTrawler
RUN R CMD INSTALL /msTrawler --no-staged-install

WORKDIR /workspace
