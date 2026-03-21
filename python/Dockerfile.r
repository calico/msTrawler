FROM python:3.12-slim

WORKDIR /app

# Install system dependencies for Python + R + rpy2
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    git \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    cmake \
    libnlopt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('reshape2', 'nlme', 'nloptr', 'lme4', 'car', 'pbkrtest'), repos='https://cloud.r-project.org')"

# Copy full repo
COPY . /repo

# Install msTrawler R package (needed by the R backend)
RUN R CMD INSTALL /repo --no-staged-install

# Install Python package with R backend
WORKDIR /repo/python
RUN pip install --no-cache-dir ".[r,dev]"

WORKDIR /repo/python

# Default: run tests
CMD ["pytest", "tests/", "-v"]
