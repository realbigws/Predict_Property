
# make the subdir first
CXX = g++
FLAG = -O3
SHELL :=/bin/bash

# compile directory
all:
	$(CXX) DeepCNF_SAS_Con.cpp DeepCNF.cpp Chebyshev.cpp DeepCNF_Misc.cpp getopt.cpp profile.cpp seq.cpp $(FLAG) -o DeepCNF_SAS_Con
	mv DeepCNF_SAS_Con ../
	
# clean all the temp files
clean:
