#!/bin/bash

set -e

VERSION=${1:-"master"}
BUILD_JOBS=${2:-16}
VERSION_TAG=${3:-""}

export CMAKE_GENERATOR=Ninja

declare -A python_versions
#python_versions[cp36-cp36m]=/opt/python/cp36-cp36m/bin/python
#python_versions[cp37-cp37m]=/opt/python/cp37-cp37m/bin/python
#python_versions[cp38-cp38]=/opt/python/cp38-cp38/bin/python
#python_versions[cp39-cp39]=/opt/python/cp39-cp39/bin/python
#python_versions[cp310-cp310]=/opt/python/cp310-cp310/bin/python
#python_versions[cp311-cp311]=/opt/python/cp311-cp311/bin/python
python_versions[cp312-cp312]=/opt/python/cp312-cp312/bin/python

for python_bin in ${python_versions[*]}
do
  ${python_bin} -m pip install pip --upgrade
  ${python_bin} -m pip install wheel setuptools twine pytest-runner auditwheel scikit-build cmake numpy
done

DIR=`pwd`
# Setup opm modules
git clone https://github.com/OPM/opm-utilities
VERSION=py_gw2
git clone https://github.com/hakonhagland/opm-common -b $VERSION
git clone https://github.com/hakonhagland/opm-grid -b $VERSION
git clone https://github.com/hakonhagland/opm-simulators-1 -b $VERSION opm-simulators

ln -sf opm-utilities/opm-super/CMakeLists.txt CMakeLists.txt
sed -e 's/add_subdirectory(opm-upscaling)//' -e 's/add_dependencies(opmupscaling opmgrid)//g' -i CMakeLists.txt

mkdir -p /tmp/opm/wheelhouse
