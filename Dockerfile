FROM python:3 AS builder

RUN mkdir -p /etc/apt/keyrings
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc
RUN echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list

RUN apt-get update && apt-get install -y samtools freebayes temurin-8-jdk

COPY $PWD /bam-matcher
WORKDIR /bam-matcher

RUN wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.6.jar -O /bam-matcher/VarScan.jar

RUN wget "https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2" -O gatk.tar.bz2
RUN tar xjf gatk.tar.bz2 --no-same-owner && rm -f gatk.tar.bz2 && mv GenomeAnalysisTK-*/* . && rm -rf GenomeAnalysisTK-*

RUN pip install -r /bam-matcher/requirements.txt

ENTRYPOINT [ "/bam-matcher/bam-matcher.py" ]
