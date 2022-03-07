FROM dgarrimar/manta-sim2:v1.0.0

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.0.0" \
      description "Docker image for manta-sim"

## Install R:mvnfast
RUN R -e 'install.packages("mvnfast", repos="http://cloud.r-project.org/")'
