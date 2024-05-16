FROM python:3.6-slim
WORKDIR /MutaCLASH
COPY . .
RUN set -xe \
    && apt-get update \
    && apt-get install -y wget make samtools bowtie2 libgomp1 \
    && pip install -r requirements.txt
CMD [ "/bin/bash" ]