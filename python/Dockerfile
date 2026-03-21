FROM python:3.12-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy the full repo (needed for setuptools_scm + test data)
COPY . /repo

# Install the Python package
WORKDIR /repo/python
RUN pip install --no-cache-dir ".[dev]"

WORKDIR /repo/python

# Default: run tests
CMD ["pytest", "tests/", "-v"]
