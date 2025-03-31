# Dockerfile to generate PyPI packages. Needs to be run from the opm-common root folder
# Example use:
# sudo docker build -t manylinux2014_opm:built . -f python/Dockerfile

FROM quay.io/pypa/manylinux_2_28_x86_64 AS stage1
ARG version="master"
ARG build_jobs=4
ARG version_tag=""
WORKDIR /tmp/opm
COPY *.sh .
RUN /bin/bash setup-docker-image.sh
RUN /bin/bash build-deps.sh $build_jobs
RUN /bin/bash generate-pypi-package.sh "$version" "$build_jobs" "$version_tag"

