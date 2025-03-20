PROGRAM       = main

version       = JTKT
CXX           = g++
CXXFLAGS      = -g -O3 -Wall -Wno-deprecated -D$(version)
LD            = g++
LDFLAGS       = 

# Add JCorran library
# CXXFLAGS     += -I/projappl/project_2003154/lib/JCORRAN/include
# LDFLAGS      += -L/projappl/project_2003154/lib/JCORRAN/lib -lJCORRANbase
CXXFLAGS     += $(shell jcorran-config-inc)
LDFLAGS      += $(shell jcorran-config-lib)
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS += $(shell root-config --libs)

# HDRSDICT = AliJBaseTrack.h AliJFFlucAnalysis.h \
# 		AliAnalysisAnaTwoMultiCorrelations.h AliAnalysisSPCRun2.h \
# 		AliAnalysisPtVn.h

MAIN = main.cpp


all: $(PROGRAM)

$(PROGRAM): $(OBJS) $(MAIN)
	@echo "Linking $(PROGRAM) ..."
	$(CXX) -lPhysics -L$(PWD) $(MAIN) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)
	@ln -s $(shell jcorran-config-path)/include/jcorranDict_rdict.pcm jcorranDict_rdict.pcm
	@echo "done"

%.cxx:
%: %.cxx
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
	rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl jcorranDict_rdict.pcm

cl:
	clean $(PROGRAM)

# nanoDict.cc:
# 	@echo "Generating dictionary ..."
# 	@rm -f nanoDict.cc nanoDict.hh nanoDict.h
# 	@rootcling nanoDict.cc $(shell jcorran-config-path )/include/*.h
