#!/bin/bash
if [ $# -ne 3 ]
then
	echo "Usage: ./genACC.sh <tgt_file> <tmp_root> <home_root>"
	exit
fi

# ---- get arguments ----#
tgt_file=$1
tmp_root=$2
home_root=$3

# ---- process -----#
RaptorX_HOME=$home_root
fulnam=`basename $tgt_file`
relnam=${fulnam%.*}
ACCPred=$RaptorX_HOME/bin/AcconPred
$ACCPred $tgt_file 1 > $tmp_root/$relnam.acc

