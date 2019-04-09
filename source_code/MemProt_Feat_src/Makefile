

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	# make Diso_Feature_Make
	$(CXX) profile.cpp seq.cpp MemProt_Feat.cpp $(FLAG) -o MemProt_Feat
	$(CXX) profile.cpp seq.cpp MemProt_Feat_noprof.cpp $(FLAG) -o MemProt_Feat_noprof
	mv MemProt_Feat MemProt_Feat_noprof ../

# clean all the temp files
clean:
