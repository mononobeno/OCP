FROM debian:12-slim

ENV DEBIAN_FRONTEND=noninteractive

# 基本ツール + GROMACS + OpenBabel + ACPYPE
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      bash \
      coreutils \
      wget \
      sqlite3 \
      gromacs \
      openbabel \
      python3 \
      python3-pip \
      python3-openbabel \
      libarpack2 \
      liblapack3 \
      libblas3 \
      libhdf4-0-alt \
      libnetcdf19 \
      tzdata && \
    rm -rf /var/lib/apt/lists/*

# ACPYPE をインストール
RUN pip3 install --no-cache-dir acpype --break-system-packages

WORKDIR /app

ENV WORKSPACE=/workspace \
    DB_PATH=/workspace/catalog/db/ocp_results.sqlite \
    GMX_MAXBACKUP=-1

COPY gromacs_prep.sh /app/gromacs_prep.sh
RUN chmod +x /app/gromacs_prep.sh

ENTRYPOINT ["/app/gromacs_prep.sh"]
