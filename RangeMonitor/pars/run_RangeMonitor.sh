#! /bin/sh

source ~/opt/gds/etc/gds-user-env.sh

export ROOTSYS=/usr/
#export LIGOSMPART=LHO_Data
export DMTINPUT=/archive/frames/O1/raw/H1/H-H1_R-11317/H-H1_R-1131753*
#export DMTINPUT=/archive/frames/O1/raw/H1/H-H1_R-11293/H-H1_R-1129383*
#export DMTHTMLOUT=${PWD}/output

#mkdir -p ${DMTHTMLOUT}

RangeMonitor -config RangeMonitor_CAL_H1.conf 
