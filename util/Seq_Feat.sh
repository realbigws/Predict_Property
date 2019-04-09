#!/bin/bash

# ----- usage ------ #
usage()
{
	echo "Sequence Feature Generate "
	echo ""
	echo "USAGE:  ./Seq_Feat.sh <-i input_fasta> [-o out_root] [-H home] "
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta :  input protein sequence file in FASTA format "
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_root    :  default output would be XXXX.feat_noprof at the current directory. "
	echo "                  [default = './' ] "
	echo ""
	echo "-H home        :  home directory of Seq_Feat.sh [default = ~/GitBucket/Sheng_Property]"
	echo ""
	exit 1
}

if [ $# -lt 1 ];
then
        usage
fi
curdir="$(pwd)"

# ----- main directory ---#
util=bin

# ----- get arguments ----- #
#-> required arguments
input=""
input_fasta=""
Keep_file=0
out_root="./"   #-> output to current directory
home=~/GitBucket/Sheng_Property         #-> home directory

#-> parse arguments
while getopts ":i:o:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input=$OPTARG
		;;
	#-> optional arguments
	o)
		out_root=$OPTARG
		;;
	H)
		home=$OPTARG
		;;
	#-> default
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done

# ------ judge fasta or tgt -------- #
filename=`basename $input`
extension=${filename##*.}
filename=${filename%.*}
input_fasta=$input

# ------ check required arguments ------ #
#-> check input_fasta
has_fasta=0
if [ ! -f "$curdir/$input_fasta" ]
then
	if [ ! -f "$input_fasta" ]
	then
		has_fasta=0
	else
		has_fasta=1
	fi
else
	input_fasta=$curdir/$input_fasta
	has_fasta=1
fi

# ------ final check ------#
if [ $has_fasta -eq 0 ] 
then
	echo "input_fasta $input_fasta not found" >&2
	exit 1
fi

# ------ part 0 ------ # related path
if [ $has_fasta -eq 1 ]
then
	fulnam=`basename $input_fasta`
	relnam=${fulnam%.*}
fi


# ------ check home directory ---------- #
if [ ! -d "$home" ];
then
	echo "home directory $home not exist " >&2
	exit 1
fi
#-> change to absolute
if [ ! -d "$curdir/$home" ]
then
if [ ! -d "$home" ]
	then
		echo "home $home not found !!" >&2
		exit 1
	fi
else
	home=$curdir/$home
fi
#echo "home=$home"


# ------ check output directory -------- #
dir_out_root=`dirname $out_root`
nam_out_root=`basename $out_root`
if [ "$dir_out_root" == "." ]
then
	if [ "$nam_out_root" == "." ]
	then
		out_root=$curdir
	else
		out_root=$curdir/$nam_out_root
	fi
fi
mkdir -p $out_root
#-> change to absolute
if [ ! -d "$curdir/$out_root" ]
then
	if [ ! -d "$out_root" ]
	then
		echo "outroot $out_root not found !!" >&2
		exit 1
	fi
else
	out_root=$curdir/$out_root
fi


# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp="${out_root}/TMP_SEQFEAT_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp


# ----- pre process ------ #
rm -f $tmp/$relnam.seq
if [ $has_fasta -eq 1 ]
then
	cp $input_fasta $tmp/$relnam.seq
fi

# ----- main procedure ------ #
program_suc=1
for ((i=0;i<1;i++))
do
	# ----- generate predicted SSE and ACC ----- #
	$home/util/runxxxpred_single.sh $tmp/$relnam.seq $tmp $home 1> $tmp/$relnam.ws1 2> $tmp/$relnam.ws2
	rm -f $tmp/$relnam.ss $tmp/$relnam.horiz $tmp/$relnam.ws1 $tmp/$relnam.ws2
	# ----- generate feature ----- #
	$home/$util/Diso_Feature_Make_noprof $tmp/$relnam.seq $tmp/$relnam.ss2 $tmp/$relnam.solv -1 > $out_root/$relnam.feat_noprof
	OUT=$?
	if [ $OUT -ne 0 ]
	then
		echo "Failed in generating feature file (no_profile mode) for sequence $relnam"
		program_suc=0
		break
	fi
done

# ----- return back ----#
if [ $Keep_file -eq 0 ]
then
	rm -rf $tmp/
fi

# ---- exit ----- #
if [ $program_suc -ne 0 ]
then
	exit 0
else
	exit 1
fi

