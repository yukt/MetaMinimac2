#!/bin/sh

rm -rf cget/ release-build/ install.log
echo -e "Installing Dependencies - Libstatgen ..."
cget ignore xz
cget install -f requirements.txt &> install.log
mkdir release-build
cd release-build/
echo -e "Generating MakeFiles ..."
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
make
echo "Binary created at /release-build/MetaMinimac2"
