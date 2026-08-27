PROGRAM       = main
PROGRAM_HDF5  = main_hdf5

CXX           = g++
CXXFLAGS      = -g -O3 -Wall -Wno-deprecated

CXXFLAGS     += $(shell jcorran-config-inc)
LDFLAGS      += $(shell jcorran-config-lib)
SOFLAGS       = -shared

CXXFLAGS     += $(shell root-config --cflags)
LDFLAGS      += $(shell root-config --libs)

HDF5_CXX      ?= $(shell pkg-config --cflags hdf5 2>/dev/null)
HDF5_LIBS     ?= $(shell pkg-config --libs hdf5 2>/dev/null) -lhdf5_cpp

JCORRAN_PCM   = $(shell jcorran-config-path)/include/JCORRANdict_rdict.pcm

.PHONY: all main main_hdf5 clean cl

all: $(PROGRAM) $(PROGRAM_HDF5)

$(PROGRAM): main.cpp
	@echo "Linking $(PROGRAM) ..."
	$(CXX) -lPhysics -L$(PWD) main.cpp $(CXXFLAGS) $(LDFLAGS) -o $(PROGRAM)
	@$(MAKE) --no-print-directory link_pcm
	@echo "done"

$(PROGRAM_HDF5): main_hdf5.cpp
	@echo "Linking $(PROGRAM_HDF5) ..."
	@if [ -z "$(HDF5_LIBS)" ]; then \
		echo "HDF5 not found. Install HDF5 or set HDF5_CXX/HDF5_LIBS."; \
		exit 1; \
	fi
	$(CXX) -lPhysics -L$(PWD) main_hdf5.cpp $(CXXFLAGS) $(HDF5_CXX) $(LDFLAGS) $(HDF5_LIBS) -o $(PROGRAM_HDF5)
	@echo "done"

link_pcm:
	@[ -e JCORRANdict_rdict.pcm ] || ln -s $(JCORRAN_PCM) JCORRANdict_rdict.pcm

clean:
	rm -f *Dict* $(PROGRAM) $(PROGRAM_HDF5) *.o JCORRANdict_rdict.pcm

cl: clean
