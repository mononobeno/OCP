FROM debian:12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      sqlite3 \
      bash \
      coreutils \
      tzdata && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

ENV DB_PATH=/zinc-library/db/ocp_results.sqlite \
    WORKSPACE=/workspace \
    TARGET_CODE=aps_apoh \
    LIBRARY_CODE=zinc_2d_smi_v1 \
    NUM_LIGANDS=50

COPY ligand_selector.sh /app/ligand_selector.sh
RUN chmod +x /app/ligand_selector.sh

ENTRYPOINT ["/app/ligand_selector.sh"]
