BH_LIB_DIR = $(prefix)/lib
BH_PATH_SOURCE = $(abs_srcdir)/../src
BH_PATH_OBJ = ../src
INCLUDE =  -I$(BH_PATH_SOURCE) -I$(BH_PATH_SOURCE)/../.. -I$(BH_PATH_SOURCE)/../FeynDiagram/
MY_PROGRAMS_LIBS = $(shell qd-config --libs )



#ADDITIONAL_OBJ =  
# the variable ADDITIONAL_OBJ contain any additional objects you want to have linked with the programs you generate in the my_programs directory
#ADDITIONAL_OBJ = $(BH_PATH_OBJ)/../../../FeynDiagram/libFD.a 
#ADDITIONAL_OBJ = ampl_loop_qqggll.o ampl_loop_qqgll.o recursive_anal_ll.o 
ADDITIONAL_OBJ = -ldl
LIBSTOLINK = \
	$(BH_PATH_OBJ)/.libs/libtools.la \
	$(BH_PATH_OBJ)/Interface/.libs/libInterface.la \
	$(BH_PATH_OBJ)/.libs/libassembly.la \
	$(BH_PATH_OBJ)/.libs/libOLA.la \
	$(BH_PATH_OBJ)/.libs/libRatPart.la \
	$(BH_PATH_OBJ)/.libs/libCutPart.la \
	$(BH_PATH_OBJ)/ratext/.libs/libratext.la \
	$(BH_PATH_OBJ)/.libs/libBH.la \
	$(BH_PATH_OBJ)/trees_eval/.libs/libAmpl_eval.la \
	$(BH_PATH_OBJ)/rational_eval/.libs/libRateval.la \
	$(BH_PATH_OBJ)/cut/.libs/libCut_wCI.la \
	$(BH_PATH_OBJ)/cut_eval/.libs/libCuteval.la \
	$(BH_PATH_OBJ)/.libs/libBG.la       \
	$(BH_PATH_OBJ)/.libs/libBHcore.la    \
	$(BH_PATH_OBJ)/.libs/libIntegrals.la    \
	$(BH_PATH_OBJ)/.libs/libSpinors.la           

if ! PUBLIC
LIBSTOLINK += \
	$(BH_PATH_OBJ)/../FeynDiagram/.libs/libfd.la      
endif

if BH_USE_GMP
LIBSTOLINK += $(BH_PATH_OBJ)/.libs/libmpreal.la
endif

LINKER_OPTS = \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/Interface/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/cut/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/cut_eval/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/ratext/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/rational_eval/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/../FeynDiagram/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/trees/.libs \
	-Wl,-rpath -Wl,$(BH_PATH_OBJ)/trees_eval/.libs


if BHDEBUG
     CXXFLAGS += -DBH_DEBUG_ON
endif

if BHTIMING
     CXXFLAGS += -DBH_TIMING_ON
endif

LINKER_ARGS = $(LIBSTOLINK)  $(ADDITIONAL_OBJ)

if BH_USE_GMP
LINKER_ARGS += -lmpfr -lgmpxx -lgmp
CXXFLAGS += -DBH_USE_GMP
endif 

include $(DEPDIR)/*.Po

%.o:	$(abs_srcdir)/%.cpp 
	test -f $(DEPDIR)/$*.Po || mkdir -p .deps ; echo "#dummy" > $(DEPDIR)/$*.Po 
	$(CXX) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@  $< 
	mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

LINKER_FLAGS=$(AM_LDFLAGS) 
if GMP_PATH_SET 
   LINKER_FLAGS += -L$(AM_GMP_PATH)
endif

%:	%.o $(LIBSTOLINK)
	libtool --mode=link --tag=CXX   $(CXX) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -o $* $*.o $(LINKER_FLAGS) $(LINKER_OPTS) $(LINKER_ARGS)

.SECONDARY: $(ADDITIONAL_OBJ) 
