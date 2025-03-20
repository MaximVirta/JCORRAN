PROGRAM       = main

CXX           = g++
CXXFLAGS      = -g -O3 -Wall -Wno-deprecated

# Add JCorran library
CXXFLAGS     += $(shell jcorran-config-inc)
LDFLAGS      += $(shell jcorran-config-lib)
SOFLAGS       = -shared

# Add ROOT dependencies
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS += $(shell root-config --libs)

MAIN = main.cpp

all: $(PROGRAM)

# Link the program and jcorranDict file
$(PROGRAM): $(MAIN)
	@echo "Linking $(PROGRAM) ..."
	$(CXX) -lPhysics -L$(PWD) $(MAIN) $(CXXFLAGS) $(LDFLAGS) -o $(PROGRAM)
	@ln -s $(shell jcorran-config-path)/include/jcorranDict_rdict.pcm jcorranDict_rdict.pcm # IMPORTANT
	@echo "done"

# Removal of files 
clean:
	rm -f *Dict* $(PROGRAM).o $(PROGRAM)

cl:
	clean $(PROGRAM)
