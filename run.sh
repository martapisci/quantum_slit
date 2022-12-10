#!/bin/bash
M=$1
DT=$2
T=$3
XC=$4
YC=$5
SIGMAX=$6
SIGMAY=$7
PX=$8
PY=$9
V0=${10}
NSLIT=${11}

sed -i.bak "7s/.*/$M/g;10s/.*/$DT/g;13s/.*/$T/g;16s/.*/$XC/g;19s/.*/$YC/g;22s/.*/$SIGMAX/g;25s/.*/$SIGMAY/g;28s/.*/$PX/g;31s/.*/$PY/g;34s/.*/$V0/g;37s/.*/$NSLIT/g;" "src/parameters.txt"

cd build
make
./probability
cd ..