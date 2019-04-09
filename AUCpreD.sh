#!/bin/bash


# ----- usage ------ #
usage()
{
	echo "AUCpreD v1.04 [Oct-01-2018] "
	echo "    Predict order/disorder regions using sequence or profile information"
	echo ""
	echo "USAGE:  ./AUCpreD.sh <-i input_fasta | input_tgt> [-o out_root] "
	echo "                     [-t threshold] [-k keep_file] [-l real_label] [-H home] "
	echo ""
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta :  input protein sequence file in FASTA format"
	echo "(or)"
	echo "-i input_tgt   :  input protein profile file in TGT format"
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_root    :  default output would be XXXX.diso_MODE at the current directory,"
	echo "                  where XXXX is the input name, and MODE is profile or noprof."
	echo "                  [default = './' ]"
	echo ""
	echo "-t threshold   :  threshold to determine disordered residue. [default = 0.5]"
	echo ""
	echo "-k keep_file   :  keep the intermediate files if its value is 1 [default = 0]"
	echo ""
	echo "-l real_label  :  real Order/Disorder label file, in three lines [default = null]"
	echo ""
	echo "-H home        :  home directory of AUCpreD.sh "
	echo "                  [default = `dirname $0`]"
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
input_tgt=""
amino_only=0
label_file=""
out_root="./"   #-> output to current directory
#-> optional arguments
threshold=""
threshold_ami=0.5
threshold_pro=0.5
Keep_file=0
home=`dirname $0`         #-> home directory


#-> parse arguments
while getopts ":i:o:t:k:l:H:" opt;
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
	t)
		threshold=$OPTARG
		;;
	k)
		Keep_file=$OPTARG
		;;
	l)
		label_file=$OPTARG
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
if [ "$extension" == "tgt" ]
then
	input_tgt=$input
	amino_only=0
else
	input_fasta=$input
	amino_only=1
fi

# ------ check required arguments ------ #
#-> check input_tgt
has_tgt=0
if [ ! -f "$curdir/$input_tgt" ]
then
	if [ ! -f "$input_tgt" ]
	then
		has_tgt=0
	else
		has_tgt=1
	fi
else
	input_tgt=$curdir/$input_tgt
	has_tgt=1
fi

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
if [ $has_fasta -eq 0 ] && [ $has_tgt -eq 0 ]
then
	echo "input_fasta $input_fasta or input_tgt $input_tgt not found" >&2
	exit 1
fi

# ------ part 0 ------ # related path
if [ $has_fasta -eq 1 ]
then
	fulnam=`basename $input_fasta`
	relnam=${fulnam%.*}
fi
if [ $has_tgt -eq 1 ]
then
	fulnam=`basename $input_tgt`
	relnam=${fulnam%.*}
fi

# ------ check home directory ---------- #
# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check output directory ------#
out_root=`readlink -f $out_root`
mkdir -p $out_root
out_root=`readlink -f $out_root`

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp="${out_root}/TMP_AUCPRED_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp

# ----- pre process ------ #
rm -f $tmp/$relnam.seq
if [ $has_fasta -eq 1 ]
then
	cp $input_fasta $tmp/$relnam.seq
fi
if [ $has_tgt -eq 1 ]
then
	cp $input_tgt $tmp/$relnam.tgt
	echo ">$relnam" > $tmp/$relnam.seq
	head -n4 $tmp/$relnam.tgt | tail -n1 | awk '{print $3}' >> $tmp/$relnam.seq
fi

# ----- main procedure ------ #
program_suc=1
for ((i=0;i<1;i++))
do
	# ------------ profile mode ---------- #
	if [ $amino_only -eq 0 ]
	then
		# ----- buildFeature ------ #
		if [ ! -f "$tmp/$relnam.tgt" ]
		then
			echo "please generate TGT file first by TGT_Package for $relnam"
			program_suc=0
			break
		fi
		# ----- DeepCNF_SS_Con ----- #
		$home/$util/DeepCNF_SS_Con -t $tmp/$relnam.tgt > $tmp/$relnam.ss8
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in generating SS8 file for sequence $relnam"
			program_suc=0
			break
		fi
		# ----- AcconPred ----- #
		$home/$util/AcconPred $tmp/$relnam.tgt 1 > $tmp/$relnam.acc
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in generating ACC file for sequence $relnam"
			program_suc=0
			break
		fi
		# ----- generate feature ----- #
		if [ "$label_file" == "" ]
		then
			$home/$util/Diso_Feature_Make $tmp/$relnam.tgt $tmp/$relnam.ss8 $tmp/$relnam.acc -1 > $tmp/$relnam.feat_profile
		else
			$home/$util/Diso_Feature_Make $tmp/$relnam.tgt $tmp/$relnam.ss8 $tmp/$relnam.acc $label_file > $tmp/$relnam.feat_profile
		fi
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in generating feature file (profile mode) for sequence $relnam"
			program_suc=0
			break
		fi
		# ----- determine threshold --- #
		if [ "$threshold" == "" ]
		then
			threshold=$threshold_pro
		fi

	# ------------ amionly mode ---------- #
	else
		# ----- generate predicted SSE and ACC ----- #
		$home/util/runxxxpred_single.sh $tmp/$relnam.seq $tmp $home 1> $tmp/$relnam.ws1 2> $tmp/$relnam.ws2
		rm -f $tmp/$relnam.ss $tmp/$relnam.horiz $tmp/$relnam.ws1 $tmp/$relnam.ws2
		# ----- generate feature ----- #
		if [ "$label_file" == "" ]
		then
			$home/$util/Diso_Feature_Make_noprof $tmp/$relnam.seq $tmp/$relnam.ss2 $tmp/$relnam.solv -1 > $tmp/$relnam.feat_noprof
		else
			$home/$util/Diso_Feature_Make_noprof $tmp/$relnam.seq $tmp/$relnam.ss2 $tmp/$relnam.solv $label_file > $tmp/$relnam.feat_noprof
		fi
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in generating feature file (no_profile mode) for sequence $relnam"
			program_suc=0
			break
		fi
		# ----- determine threshold --- #
		if [ "$threshold" == "" ]
		then
			threshold=$threshold_ami
		fi
	fi

	# ---------- predict order/disorder regions ----------- #
	outnam=$relnam.diso
	if [ $amino_only -eq 0 ]
	then
		$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_profile -w 5,5 -d 50,50 -s 3 -l 148 -m $home/parameters/AUCpreD_profile_model \
			> $tmp/$relnam.diso_profile 2> $tmp/$relnam.pred_log2
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in prediction of order/disorder (profile mode )for sequence $relnam"
			program_suc=0
			break
		fi
		outnam=$relnam.diso_profile
	else
		$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_noprof -w 5,5 -d 50,50 -s 3 -l 87 -m $home/parameters/AUCpreD_noprof_model \
			> $tmp/$relnam.diso_noprof 2> $tmp/$relnam.pred_log2
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in prediction of order/disorder (no_profile mode) for sequence $relnam"
			program_suc=0
			break
		fi
		outnam=$relnam.diso_noprof
	fi

	# ------ final copy ------ #
	cp $tmp/$outnam $out_root/$relnam.diso_prev


	# ---------- transform the raw result into DisoPred format ------ #
	$home/$util/DisoPred_Trans $tmp/$relnam.seq $tmp/$outnam $threshold $amino_only > $out_root/$outnam

done

# ----- return back ----#
if [ $Keep_file -eq 0 ]
then
	rm -rf $tmp/
else
	rm -rf $out_root/"TMP_AUCPRED_"${relnam}
	mv $tmp $out_root/"TMP_AUCPRED_"${relnam}
fi

# ---- exit ----- #
if [ $program_suc -ne 0 ]
then
	exit 0
else
	exit 1
fi

