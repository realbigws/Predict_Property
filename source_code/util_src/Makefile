

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:

	# ---- for TGT_Upgrede ----- #
	$(CXX) TPL_TGT_Get_PSM.cpp $(FLAG) -o TPL_TGT_Get_PSM
	$(CXX) TPL_TGT_Get_PSP.cpp $(FLAG) -o TPL_TGT_Get_PSP	
	$(CXX) TPL_TGT_Get_HMM.cpp $(FLAG) -o TPL_TGT_Get_HMM
	$(CXX) TPL_Get_SS2.cpp $(FLAG) -o TPL_Get_SS2
	$(CXX) TGT_Get_SS2.cpp $(FLAG) -o TGT_Get_SS2
	$(CXX) TPL_Get_SEQ.cpp $(FLAG) -o TPL_Get_SEQ
	$(CXX) TGT_Get_SEQ.cpp $(FLAG) -o TGT_Get_SEQ
	# move to ../ directory
	mv TPL_TGT_Get_PSM TPL_TGT_Get_PSP TPL_TGT_Get_HMM TPL_Get_SS2 TGT_Get_SS2 TPL_Get_SEQ TGT_Get_SEQ ../

	# ------ for RaptorX_Property -------#
	# Verify_FASTA
	$(CXX) Verify_FASTA.cpp $(FLAG) -o Verify_FASTA
	# DisoPred_Trans
	$(CXX) DisoPred_Trans.cpp $(FLAG) -o DisoPred_Trans
	# TopoPred_Trans
	$(CXX) TopoPred_Trans.cpp $(FLAG) -o TopoPred_Trans
	# SS3/SS8/ACC Trans
	$(CXX) Label_Parser.cpp $(FLAG) -o Label_Parser
	# move to ../ directory
	mv Verify_FASTA DisoPred_Trans TopoPred_Trans Label_Parser ../ 

# clean all the temp files
clean:
