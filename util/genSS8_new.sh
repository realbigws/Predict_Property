#!/bin/bash
if [ $# -ne 3 ]
then
	echo "Usage: ./genSS8.sh <tgt_file> <tmp_root> <home_root>"
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
SS8Pred=$RaptorX_HOME/bin/DeepCNF_SS_Con
$SS8Pred -t $tgt_file > $tmp_root/$relnam.ss8

