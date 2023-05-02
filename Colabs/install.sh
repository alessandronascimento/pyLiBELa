#!/bin/bash
apt-get install libboost-python1.71-dev openbabel libopenbabel-dev python-dev-is-python3 zlib1g-dev libeigen3-dev libgsl-dev libnlopt-cxx-dev libgsl-dev
make clean
make FILES
make -j4
