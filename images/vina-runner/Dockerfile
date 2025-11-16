FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      wget \
      sqlite3 \
      jq \
      autodock-vina && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY vina_runner.sh /app/vina_runner.sh

RUN chmod +x /app/vina_runner.sh

# デフォルトのエントリポイント
ENTRYPOINT ["/app/vina_runner.sh"]
