

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	# make Diso_Feature_Make
	$(CXX) profile.cpp seq.cpp Diso_Feature_Make.cpp $(FLAG) -o Diso_Feature_Make
	$(CXX) profile.cpp seq.cpp Diso_Feature_Make_noprof.cpp $(FLAG) -o Diso_Feature_Make_noprof
	mv Diso_Feature_Make Diso_Feature_Make_noprof ../

# clean all the temp files
clean:
