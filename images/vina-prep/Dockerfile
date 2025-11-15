FROM debian:12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      bash \
      coreutils \
      sqlite3 \
      openbabel \
      tzdata && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

ENV WORKSPACE=/workspace \
    TARGET_CODE=aps_apoh \
    LIBRARY_CODE=zinc_2d_smi_v1 \
    DB_PATH=/zinc-library/db/ocp_results.sqlite

COPY vina_prep.sh /app/vina_prep.sh
RUN chmod +x /app/vina_prep.sh

ENTRYPOINT ["/app/vina_prep.sh"]
