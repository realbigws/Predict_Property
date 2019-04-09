#!/bin/bash

# This is a simple script which will carry out all of the basic steps
# required to make a PSIPRED prediction. Note that it assumes that the
# following programs are available in the appropriate directories:
# seq2mtx - PSIPRED V3 program
# psipred - PSIPRED V3 program
# psipass2 - PSIPRED V3 program

# NOTE: Script modified to be more cluster friendly (DTJ April 2008)

if [ $# -ne 3 ]
then
	echo "Usage: ./runxxxpred_single.sh <input_mtx> <out_dir> <home_root> "
	exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

INPUTMTX=$1
DESTDIR=$2
home_root=$3
PSIPREDDIR=$home_root/util/psisolvpred


###### BY TINA, May 6, 2003
# first make PSP directory if it doesn't exist
if [ ! -d $DESTDIR ] ; then
    mkdir $DESTDIR
fi
###### BY TINA, May 6, 2003


# Where the PSIPRED V3 programs have been installed
execdir=$PSIPREDDIR/bin
# Where the PSIPRED V3 data files have been installed
datadir=$PSIPREDDIR/data


#----- we must create a unique name here !!!! ------- ## 2018.10.22 by Sheng Wang
fulnam=`basename $1`
bname=${fulnam%.*}
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
rootname="TMP_PSIPRED_${RANDOM}_${DATE}_${bname}"
#----- we must create a unique name here !!!! ------- ## over



# ----- run psipred ------ #
echo "Generating mtx file from sequence" $1 "..."
$execdir/seq2mtx $1 > $DESTDIR/$rootname.mtx
if [ $? -ne 0 ] 
then
    echo "FATAL: Error whilst running makemat - script terminated!"
    exit 1
fi

echo "Predicting secondary structure based on single sequence ..."
echo Pass1 ...
$execdir/psipred $DESTDIR/$rootname.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $DESTDIR/$rootname.ss
if [ $? -ne 0 ]
then
    echo "FATAL: Error whilst running psipred - script terminated!"
    exit 1
fi

echo Pass2 ...
$execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.0 $DESTDIR/$rootname.ss2 $DESTDIR/$rootname.ss > $DESTDIR/$rootname.horiz
if [ $? -ne 0 ]
then
    echo "FATAL: Error whilst running psipass2 - script terminated!"
    exit 1
fi


# Run SOLVPRED
$execdir/solvpred $DESTDIR/$rootname.mtx $datadir/weights_solv.dat > $DESTDIR/$rootname.solv
if [ $? -ne 0 ]
then
    echo "FATAL: Error whilst running MetaPSICOV solvpred - script terminated!"
    exit 1
fi


echo "Final output files:" $rootname.ss2 $rootname.horiz $rootname.ss $rootname.solv
mv $DESTDIR/$rootname.ss $DESTDIR/$bname.ss
mv $DESTDIR/$rootname.ss2 $DESTDIR/$bname.ss2
mv $DESTDIR/$rootname.horiz $DESTDIR/$bname.horiz
mv $DESTDIR/$rootname.solv $DESTDIR/$bname.solv


# Remove temporary files
echo Cleaning up ...
rm -f $DESTDIR/$rootname.mtx

