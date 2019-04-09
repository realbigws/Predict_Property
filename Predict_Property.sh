#!/bin/bash


# ----- usage ------ #
usage()
{
	echo "Predict_Property v1.02 [Mar-05-2019] "
	echo "    Predict protein local properties using sequence or profile information"
	echo ""
	echo "USAGE:  ./Predict_Property.sh <-i input_fasta | input_tgt> [-o out_root]"
	echo "                              [-t diso_thres] [-T topo_thres] [-H home] "
	echo ""
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta :  input protein sequence file in FASTA format"
	echo "(or)"
	echo "-i input_tgt   :  input protein profile file in TGT format"
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_root    :  default output would the current directory. "
	echo "                  [default = './\${input_name}_PROP'] "
	echo ""
	echo "-t diso_thres  :  threshold to determine disordered residue. [default = 0.5]"
	echo ""
	echo "-T topo_thres  :  threshold to determine transmembrane residue. [default = 0.5]"
	echo ""
	echo "-H home        :  home directory of Predict_Property.sh "
	echo "                  [default = `dirname $0`] "
	echo ""
	exit 1
}

if [ $# -lt 1 ];
then
        usage
fi
curdir="$(pwd)"


# ----- get arguments ----- #
#-> required arguments
input=""
input_fasta=""
input_tgt=""
PROF_or_NOT=1
#-> optional arguments
out_root=""   #-> output to current directory
diso_thres=0.5
topo_thres=0.5
home=`dirname $0`         #-> home directory

#-> parse arguments
while getopts ":i:o:t:T:H:" opt;
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
		diso_thres=$OPTARG
		;;
	T)
		topo_thres=$OPTARG
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



# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check input fasta ------#
if [ ! -s "$input" ]
then
	echo "input $input not found !!" >&2
	exit 1
fi
input=`readlink -f $input`
fulnam=`basename $input`
relnam=${fulnam%.*}

# ------ judge fasta or tgt -------- #
filename=`basename $input`
extension=${filename##*.}
filename=${filename%.*}
if [ "$extension" == "tgt" ]
then
	input_tgt=$input
	PROF_or_NOT=1
else
	input_fasta=$input
	PROF_or_NOT=0
fi

# ------ check output directory -------- #
if [ "$out_root" == "" ]
then
	out_root=${relnam}_PROP
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp="${out_root}/TMP_TGTGEN_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp




# ---- running ---------#
if true
then

	#------- start -------#
	util=bin
	program_suc=1
	for ((i=0;i<1;i++))
	do
		if [ $PROF_or_NOT -eq 1 ] #-> use profile
		then

			cp $input_tgt $tmp/$relnam.tgt
			$home/util/TGT_Get_SEQ $tmp/$relnam.tgt $tmp/$relnam.seq
			#-> 1.1 TGT Update
			echo "Running TGT_Update to upgrade TGT file for sequence $relnam"
			tmptmp=$out_root/TMPTMP"_"$relnam"_"$RANDOM
			mkdir -p $tmptmp
			mkdir -p $tmp/update/
			$home/$util/TGT_Update -i $tmp/$relnam.tgt -o $tmp/update/$relnam.tgt -t $tmptmp -H $home
			rm -rf $tmptmp
			#-> 2. generate SS3/SS8 file
			#--> SS8
			$home/$util/DeepCNF_SS_Con -t $tmp/$relnam.tgt -s 0 > $tmp/$relnam.ss8
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating SS8 file for sequence $relnam"
				program_suc=0
				break
			fi
			#--> SS3
			$home/$util/DeepCNF_SS_Con -t $tmp/$relnam.tgt -s 1 > $tmp/$relnam.ss3
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating SS3 file for sequence $relnam"
				program_suc=0
				break
			fi
			#-> 3. generate ACC/CN file
			#--> ACC
			$home/$util/DeepCNF_SAS_Con -t $tmp/update/$relnam.tgt -m 0 > $tmp/$relnam.acc
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating ACC file for sequence $relnam"
				program_suc=0
				break
			fi
			#--> CN
			$home/$util/AcconPred $tmp/$relnam.tgt 0 > $tmp/$relnam.cn
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating CN file for sequence $relnam"
				program_suc=0
				break
			fi
			#-> 4. generate DISO file
			$home/AUCpreD.sh -i $tmp/$relnam.tgt -o $tmp -t $diso_thres -H $home
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating DISO file for sequence $relnam"
				program_suc=0
				break
			fi
			mv $tmp/$relnam.diso_profile $tmp/$relnam.diso
			#-> 5. generate TM file
			$home/PDBTM_Topology_Pred.sh -i $tmp/$relnam.tgt -o $tmp -t $topo_thres -H $home
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating TM file for sequence $relnam"
				program_suc=0
				break
			fi
			mv $tmp/$relnam.tm2_profile $tmp/$relnam.tm2
			mv $tmp/$relnam.tm8_profile $tmp/$relnam.tm8

		else     #-> not use profile

			#-> 1. generate feature file
			$home/util/Verify_FASTA $input_fasta $tmp/$relnam.seq
			cp $input_fasta $tmp/$relnam.fasta_raw
			#--> feat_file
			$home/util/Seq_Feat.sh -i $tmp/$relnam.seq -o $tmp -H $home
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating noprof_feat file for sequence $relnam"
				program_suc=0
				break
			fi
			#--> pred_file for SS8/SS3
			$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_noprof -w 5,5,5,5,5 -d 100,100,100,100,100 -s 8 -l 87 -m $home/parameters/ss8_noprof_model > $tmp/$relnam.ss8_noprof 2> $tmp/$relnam.noprf_log3
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in prediction of SS8/SS3 (no_profile mode) for sequence $relnam"
				program_suc=0
				break
			fi
			#--> pred_file for ACC
			$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_noprof -w 5,5,5,5,5 -d 100,100,100,100,100 -s 3 -l 87 -m $home/parameters/acc_noprof_model > $tmp/$relnam.acc_noprof 2> $tmp/$relnam.noprf_log4
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in prediction of ACC (no_profile mode) for sequence $relnam"
				program_suc=0
				break
			fi
			#--> pred_file for CN
			$home/$util/DeepCNF_Pred -i $tmp/$relnam.feat_noprof -w 5,5,5,5,5 -d 100,100,100,100,100 -s 15 -l 87 -m $home/parameters/cn_noprof_model > $tmp/$relnam.cn_noprof 2> $tmp/$relnam.noprf_log5
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in prediction of CN (no_profile mode) for sequence $relnam"
				program_suc=0
				break
			fi
			#-> 2. generate SS3/SS8 file
			#--> SS8
			$home/$util/Label_Parser $tmp/$relnam.seq $tmp/$relnam.ss8_noprof 0 > $tmp/$relnam.ss8_
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating SS8 noprof_file for sequence $relnam"
				program_suc=0
				break
			fi
			#--> SS3
			$home/$util/Label_Parser $tmp/$relnam.seq $tmp/$relnam.ss8_noprof 1 > $tmp/$relnam.ss3_
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating SS3 noprof_file for sequence $relnam"
				program_suc=0
				break
			fi
			#-> 3. generate ACC/CN file
			#--> ACC
			$home/$util/Label_Parser $tmp/$relnam.seq $tmp/$relnam.acc_noprof 2 > $tmp/$relnam.acc_
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating ACC noprof_file for sequence $relnam"
				program_suc=0
				break
			fi
			#--> CN
			$home/$util/Label_Parser $tmp/$relnam.seq $tmp/$relnam.cn_noprof 3 > $tmp/$relnam.cn_
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating CN noprof_file for sequence $relnam"
				program_suc=0
				break
			fi
			#-> 4. generate DISO file
			$home/AUCpreD.sh -i $tmp/$relnam.seq -o $tmp -t $diso_thres -H $home
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating DISO noprof_file for sequence $relnam"
				program_suc=0
				break
			fi
			mv $tmp/$relnam.diso_noprof $tmp/$relnam.diso_
			#-> 5.generate TM file
			$home/PDBTM_Topology_Pred.sh -i $tmp/$relnam.seq -o $tmp -t $topo_thres -H $home
			OUT=$?
			if [ $OUT -ne 0 ]
			then
				echo "Failed in generating TM noprof_file for sequence $relnam"
				program_suc=0
				break
			fi
			mv $tmp/$relnam.tm2_noprof $tmp/$relnam.tm2_
			mv $tmp/$relnam.tm8_noprof $tmp/$relnam.tm8_

		fi
	done
	# ----------- end ------------- #
	if [ $program_suc -ne 1 ]
	then
        	exit 1
	fi
	# ----------- copy to $out_root/ ----- #
	if [ $PROF_or_NOT -eq 1 ] #-> use profile
	then

		cp $home/util/0README $out_root/0README.txt
		cp $tmp/$relnam.seq $out_root/$relnam.fasta.txt
		cp $tmp/$relnam.seq $out_root/$relnam.seq.txt
		# copy A3M and TGT if available
		if [ -f $tmp/$relnam.tgt ] 
		then 
			mkdir -p $out_root/Profile_data
			cp $tmp/$relnam.tgt $out_root/Profile_data/
		fi
		if [ -f $tmp/$relnam.a3m ]
		then
			mkdir -p $out_root/Profile_data
			cp $tmp/$relnam.a3m $out_root/Profile_data/
		fi
		# make simple prediction
		#-> all types
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			awk '{if(NF>0){print $0}}' $tmp/$relnam.${ttt} > $out_root/$relnam.${ttt}.txt
			echo ">$relnam" > $out_root/$relnam.${ttt}_simp.txt
			grep -v "#" $out_root/$relnam.${ttt}.txt | awk '{printf $2}END{printf "\n"}' >> $out_root/$relnam.${ttt}_simp.txt
			grep -v "#" $out_root/$relnam.${ttt}.txt | awk '{printf $3}END{printf "\n"}' >> $out_root/$relnam.${ttt}_simp.txt
		done
		# make overall prediction
		head -n1 $out_root/$relnam.fasta.txt > $out_root/$relnam.all.txt
		tail -n1 $out_root/$relnam.seq.txt >> $out_root/$relnam.all.txt
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			tail -n1 $out_root/$relnam.${ttt}_simp.txt >> $out_root/$relnam.all.txt
		done
		printf "\n\n" >> $out_root/$relnam.all.txt
		printf "\n\n#---------------- details of SS3 prediction ---------------------------\n" > $out_root/$relnam.all.ss3
		printf "\n\n#---------------- details of SS8 prediction ---------------------------\n" > $out_root/$relnam.all.ss8
		printf "\n\n#---------------- details of ACC prediction ---------------------------\n" > $out_root/$relnam.all.acc
		printf "\n\n#---------------- details of DISO prediction --------------------------\n" > $out_root/$relnam.all.diso
		printf "\n\n#---------------- details of TM2 prediction ---------------------------\n" > $out_root/$relnam.all.tm2
		printf "\n\n#---------------- details of TM8 prediction ---------------------------\n" > $out_root/$relnam.all.tm8
		cat $out_root/$relnam.all.txt > $out_root/$relnam.all.txt_
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			cat $out_root/$relnam.all.${ttt} $out_root/$relnam.${ttt}.txt >> $out_root/$relnam.all.txt_
		done
		mv $out_root/$relnam.all.txt_ $out_root/$relnam.all.txt
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			rm -f $out_root/$relnam.all.${ttt}
		done

	else                      #-> not use profile

		cp $home/util/0README_noprof $out_root/0README.txt
		cp $tmp/$relnam.fasta_raw $out_root/$relnam.fasta.txt
		cp $tmp/$relnam.seq $out_root/$relnam.seq.txt
		# make simple prediction
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			awk '{if(NF>0){print $0}}' $tmp/$relnam.${ttt}_ > $out_root/$relnam.${ttt}_noprof.txt
			echo ">$relnam" > $out_root/$relnam.${ttt}_noprof_simp.txt
			grep -v "#" $out_root/$relnam.${ttt}_noprof.txt | awk '{printf $2}END{printf "\n"}' >> $out_root/$relnam.${ttt}_noprof_simp.txt
			grep -v "#" $out_root/$relnam.${ttt}_noprof.txt | awk '{printf $3}END{printf "\n"}' >> $out_root/$relnam.${ttt}_noprof_simp.txt
		done
		# make overall prediction
		head -n1 $out_root/$relnam.fasta.txt > $out_root/$relnam.all.txt
		tail -n1 $out_root/$relnam.seq.txt >> $out_root/$relnam.all.txt
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			tail -n1 $out_root/$relnam.${ttt}_noprof_simp.txt >> $out_root/$relnam.all.txt
		done
		printf "\n\n" >> $out_root/$relnam.all.txt
		printf "\n\n#---------------- details of SS3 prediction ---------------------------\n" > $out_root/$relnam.all.ss3
		printf "\n\n#---------------- details of SS8 prediction ---------------------------\n" > $out_root/$relnam.all.ss8
		printf "\n\n#---------------- details of ACC prediction ---------------------------\n" > $out_root/$relnam.all.acc
		printf "\n\n#---------------- details of DISO prediction --------------------------\n" > $out_root/$relnam.all.diso
		printf "\n\n#---------------- details of TM2 prediction ---------------------------\n" > $out_root/$relnam.all.tm2
		printf "\n\n#---------------- details of TM8 prediction ---------------------------\n" > $out_root/$relnam.all.tm8
		cat $out_root/$relnam.all.txt > $out_root/$relnam.all.txt_
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			cat $out_root/$relnam.all.${ttt} $out_root/$relnam.${ttt}_noprof.txt >> $out_root/$relnam.all.txt_
		done
		mv $out_root/$relnam.all.txt_ $out_root/$relnam.all.txt
		mv $out_root/$relnam.all.txt $out_root/$relnam.all_noprof.txt
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			rm -f $out_root/$relnam.all.${ttt}
		done

	fi

	# --------- rename not use profile mode ----- #
	if [ $PROF_or_NOT -ne 1 ]
	then
		mv $out_root/0README.txt $out_root/0README_noprof
		mv $out_root/$relnam.fasta.txt $out_root/$relnam.fasta
		mv $out_root/$relnam.seq.txt $out_root/$relnam.seq
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			#--> raw prediction result
			mv $out_root/$relnam.${ttt}_noprof.txt $out_root/$relnam.${ttt}
			#--> simp prediction result
			mv $out_root/$relnam.${ttt}_noprof_simp.txt $out_root/$relnam.${ttt}_simp
		done
		#-> overall prediction result
		mv $out_root/$relnam.all_noprof.txt $out_root/$relnam.all
	else
		mv $out_root/0README.txt $out_root/0README
		mv $out_root/$relnam.fasta.txt $out_root/$relnam.fasta
		mv $out_root/$relnam.seq.txt $out_root/$relnam.seq
		for ttt in ss3 ss8 acc diso tm2 tm8
		do
			#--> raw prediction result
			mv $out_root/$relnam.${ttt}.txt $out_root/$relnam.${ttt}
			#--> simp prediction result
			mv $out_root/$relnam.${ttt}_simp.txt $out_root/$relnam.${ttt}_simp
		done
		#-> overall prediction result
		mv $out_root/$relnam.all.txt $out_root/$relnam.all
		#-> move tgt files
		cp $tmp/$relnam.tgt $out_root/$relnam.tgt
		cp $tmp/update/$relnam.tgt $out_root/$relnam.tgt2
		if [ -f $tmp/$relnam.a3m ]
		then
			cp $tmp/$relnam.a3m $out_root/$relnam.a3m
		fi
	fi
	# -------- prediction summary -------- #
	$home/$util/generate_simp_summary_file $out_root/$relnam.diso $out_root/$relnam.tm2 $out_root/$relnam.ss3 $out_root/$relnam.acc \
		$diso_thres $topo_thres $out_root/$relnam.summary
	

	# ------ remove temporary folder ----- #
	rm -rf $tmp/
fi

