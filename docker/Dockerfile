FROM dgarrimar/manta-sim3:v1.0.0

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.0.0" \
      description "Docker image for manta-sim (sim2)"

## Install R:mvnfast
RUN R -e 'install.packages("mvnfast", repos="http://cloud.r-project.org/")'
