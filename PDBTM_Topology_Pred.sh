#!/bin/bash


# ----- usage ------ #
usage()
{
	echo "PDBTM_Topology_Pred v1.05 [Mar-05-2019] "
	echo "    Predict PDBTM Topology labels using sequence or profile information "
	echo ""
	echo "USAGE:  ./PDBTM_Topology_Pred.sh <-i input_fasta | input_tgt> [-o out_root] "
	echo "                  [-t threshold] [-k keep_file] [-l real_label] [-H home]"
	echo ""
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta :  input protein sequence file in FASTA format"
	echo "(or)"
	echo "-i input_tgt   :  input protein profile file in TGT format"
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_root    :  default output would be XXXX.topo_MODE at the current directory,"
	echo "                  where XXXX is the input name, and MODE is profile or noprof."
	echo "                  we also output XXXX.topo_simp for simplified result."
	echo "                  [default = './' ]"
	echo ""
	echo "-t threshold   :  threshold to determine transmembrane residue. [default = 0.5]"
	echo ""
	echo "-k keep_file   :  keep the intermediate files if its value is 1 [default = 0]"
	echo ""
	echo "-l real_label  :  real PDBTM Topology label file, in  three lines [default = null]" 
	echo ""
	echo "-H home        :  home directory of PDBTM_Topology_Pred.sh "
	echo "                  [default = `dirname $0`] "
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
home=`dirname $0`        #-> home directory


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
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check output directory -------- #
mkdir -p $out_root
out_root=`readlink -f $out_root`

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp="${out_root}/TMP_PDBTM_${relnam}_${RANDOM}_${DATE}"
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
		# ----- generate feature ----- #
		if [ "$label_file" == "" ]
		then
			$home/$util/MemProt_Feat $tmp/$relnam.tgt $tmp/$relnam.ss8 null > $tmp/$relnam.feat_profile
		else
			$home/$util/MemProt_Feat $tmp/$relnam.tgt $tmp/$relnam.ss8 $label_file > $tmp/$relnam.feat_profile
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
			$home/$util/MemProt_Feat_noprof $tmp/$relnam.seq $tmp/$relnam.ss2 $tmp/$relnam.solv null > $tmp/$relnam.feat_noprof
		else
			$home/$util/MemProt_Feat_noprof $tmp/$relnam.seq $tmp/$relnam.ss2 $tmp/$relnam.solv $label_file > $tmp/$relnam.feat_noprof
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
	outnam=$relnam.topo
	outfeat=$relnam.feat
	if [ $amino_only -eq 0 ]
	then
		$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_profile -w 5,5,5,5,5 -d 100,100,100,100,100 -s 9 -l 68 -m $home/parameters/MemProt_profile_model > $tmp/$relnam.topo_profile 2> $tmp/$relnam.pred_log2
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in prediction of membrane topology (profile mode )for sequence $relnam"
			program_suc=0
			break
		fi
		outnam=$relnam.topo_profile
		outfeat=$relnam.feat_profile
	else
		$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_noprof -w 5,5,5,5,5 -d 100,100,100,100,100 -s 9 -l 87 -m $home/parameters/MemProt_noprof_model > $tmp/$relnam.topo_noprof 2> $tmp/$relnam.pred_log2
		OUT=$?
		if [ $OUT -ne 0 ]
		then
			echo "Failed in prediction of membrane topology (no_profile mode) for sequence $relnam"
			program_suc=0
			break
		fi
		outnam=$relnam.topo_noprof
		outfeat=$relnam.feat_noprof
	fi

	# ------ final copy ------ #
	cp $tmp/$outnam $out_root

	# ------ transform the raw result into DisoPred format ------ #
	if [ $amino_only -eq 0 ]
	then
		$home/$util/TopoPred_Trans $tmp/$relnam.seq $tmp/$outnam $threshold $amino_only 0 > $out_root/$relnam.tm2_profile
		$home/$util/TopoPred_Trans $tmp/$relnam.seq $tmp/$outnam $threshold $amino_only 1 > $out_root/$relnam.tm8_profile
	else
		$home/$util/TopoPred_Trans $tmp/$relnam.seq $tmp/$outnam $threshold $amino_only 0 > $out_root/$relnam.tm2_noprof
		$home/$util/TopoPred_Trans $tmp/$relnam.seq $tmp/$outnam $threshold $amino_only 1 > $out_root/$relnam.tm8_noprof
	fi

done



#============== output oneline prediction in FASTA format: 0 for non-TM region and 1 for TM region ===============#
if true
then
	grep -v "#" $tmp/$outnam | awk '{a=$6;if(a>0.5){printf "1"}else{printf "0"}}END{printf "\n"}' > $tmp/$relnam.topo_simp_helix
	grep -v "#" $tmp/$outnam | awk '{b=$5;if(b>0.5){printf "1"}else{printf "0"}}END{printf "\n"}' > $tmp/$relnam.topo_simp_sheet
	#-> check for predicted transmembrane region
	str1=`tail -n1 $tmp/$relnam.topo_simp_helix`
	top1="${str1//[^1]}"
	num1=${#top1}
	str2=`tail -n1 $tmp/$relnam.topo_simp_sheet`
	top2="${str2//[^1]}"
	num2=${#top2}
	#-> determine helix_type or sheet_type
	cat $tmp/$relnam.seq > $tmp/$relnam.topo_simp
	if [ $num1 -ge $num2 ]
	then
		tail -n1 $tmp/$relnam.topo_simp_helix >> $tmp/$relnam.topo_simp
	else
		tail -n1 $tmp/$relnam.topo_simp_sheet >> $tmp/$relnam.topo_simp
	fi
	#-> final copy
	cp $tmp/$relnam.topo_simp $out_root
fi


# ----- return back ----#
if [ $Keep_file -eq 0 ]
then
	rm -rf $tmp/
else
	rm -rf $out_root/"TMP_PDBTM_"${relnam}
	mv $tmp $out_root/"TMP_PDBTM_"${relnam}
fi

# ---- exit ----- #
if [ $program_suc -ne 0 ]
then
	exit 0
else
	exit 1
fi

