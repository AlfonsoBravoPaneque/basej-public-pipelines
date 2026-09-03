# ============================================================================
# protein_cpm — container for the PROTEIN_CPM process (basej-rnaproteinqc)
# ============================================================================
# Minimal Python image with pandas/numpy/pyarrow for ADT count merging and
# CPM normalization.
#
# Build (x86, local tag):
#   docker build -t basejumper_protein_cpm_0.1 .
# ============================================================================
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

COPY ./requirements.txt /requirements.txt
RUN pip3 install --no-cache-dir -r /requirements.txt

WORKDIR /data

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser && chown appuser:appuser /data
USER appuser

CMD ["/bin/bash"]
