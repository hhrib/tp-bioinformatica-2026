FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    emboss \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install biopython

WORKDIR /app

COPY src/ ./src/
COPY data/ ./data/

CMD ["bash"]
