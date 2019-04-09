

# make the subdir first
SHELL :=/bin/bash
CXX = g++
FLAG = -O3

Util_Dir = ./util_src/
DeepCNF_SS_Con_Dir = ./DeepCNF_SS_Con_src/
DeepCNF_SAS_Con_Dir = ./DeepCNF_SAS_Con_src/
AcconPred_Dir = ./AcconPred_src/
Diso_Feature_Make_Dir = ./Diso_Feature_Make_src/
Memb_Feature_Make_Dir = ./MemProt_Feat_src/
DeepCNF_Pred_Dir = ./DeepCNF_Pred_src/


all: tgt util sse acc feat pred

tgt:
	# ---- TGT_Update ----#
	$(CXX) TGT_Update_acc_ss8_new.cpp $(FLAG) -o TGT_Update
	mv TGT_Update ../bin

util:
	# ---- util ---- #
	cd $(Util_Dir)&&make
	mv Label_Parser DisoPred_Trans TopoPred_Trans ../bin
	mv Verify_FASTA ../util
	#-> for TGT_Update
	mv TPL_TGT_Get_PSM TPL_TGT_Get_PSP TPL_TGT_Get_HMM TPL_Get_SS2 TGT_Get_SS2 TPL_Get_SEQ TGT_Get_SEQ ../util
	#-> for RaptorX_Property server
	$(CXX) generate_simp_summary_file.cpp $(FLAG) -o generate_simp_summary_file
	mv generate_simp_summary_file ../bin

sse:
	# ---- DeepCNF_SS_Con ---- #
	cd $(DeepCNF_SS_Con_Dir)&&make
	mv DeepCNF_SS_Con ../bin

acc:
	# ---- AcconPred ----- #
	cd $(AcconPred_Dir)&&make
	mv AcconPred ../bin
	cd $(DeepCNF_SAS_Con_Dir)&&make
	mv DeepCNF_SAS_Con ../bin

feat:
	# ---- Diso_Feature_Make --- #
	cd $(Diso_Feature_Make_Dir)&&make
	mv Diso_Feature_Make Diso_Feature_Make_noprof ../bin
	cd $(Memb_Feature_Make_Dir)&&make
	mv MemProt_Feat MemProt_Feat_noprof ../bin

pred:
	# ---- DeepCNF_Pred ---- #
	cd $(DeepCNF_Pred_Dir)&&make
	mv DeepCNF_Pred ../bin


# clean all the temp files
clean:

